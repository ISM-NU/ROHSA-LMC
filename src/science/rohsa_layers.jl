# Velocity-layer cubes from a ROHSA decomposition.
#
# Individual ROHSA Gaussians are not unique (neighbours in velocity trade
# column from pixel to pixel), but sums over a velocity range are stable. A
# layer cube collects every Gaussian by its centre velocity relative to a
# smooth reference velocity field, so VGT can run on each layer separately.
#
# The reference is the model centroid smoothed over many beams. It follows the
# disk rotation but not the pixel-scale switches of "the strongest Gaussian",
# which would put artificial edges into the layer cubes. Membership changes
# smoothly (logistic over `width_kms`) for the same reason. The layers add up to
# the full model cube exactly.

const ROHSA_LAYER_NAMES = (:low, :main, :high)

logistic(x) = 1 / (1 + exp(-x))

"""
    velocity_layer_weights(offset_kms; split_kms=20.0, width_kms=3.0)

Weights `(low, main, high)` of a Gaussian whose centre lies `offset_kms` from
the reference velocity. They sum to one.
"""
function velocity_layer_weights(offset_kms::Real; split_kms::Real=20.0, width_kms::Real=3.0)
    width_kms > 0 || throw(ArgumentError("width_kms must be > 0"))
    low = logistic((-split_kms - offset_kms) / width_kms)
    high = logistic((offset_kms - split_kms) / width_kms)
    return (low, 1 - low - high, high)
end

"""
    reference_centroid_channels(params; sigma_pixels=30.0)

Column-weighted Gaussian centre per pixel (channel units), smoothed with a
Gaussian of `sigma_pixels`. Pixels without emission are filled by normalized
convolution.
"""
function reference_centroid_channels(params; sigma_pixels::Real=30.0)
    n_gauss = size(params, 1) ÷ 3
    nx, ny = size(params, 2), size(params, 3)
    num = zeros(Float64, nx, ny)
    den = zeros(Float64, nx, ny)
    @inbounds for y in 1:ny, x in 1:nx, g in 1:n_gauss
        w = max(params[3g - 2, x, y], 0.0) * max(params[3g, x, y], 0.0)  # ∝ integrated intensity
        num[x, y] += w * params[3g - 1, x, y]
        den[x, y] += w
    end
    kernel = Kernel.gaussian((Float64(sigma_pixels), Float64(sigma_pixels)))
    snum = imfilter(num, kernel, Fill(0.0))
    sden = imfilter(den, kernel, Fill(0.0))
    return map((n, d) -> d > 0 ? n / d : NaN, snum, sden)
end

"""
    build_layer_cubes(params, nv, channel_kms; reference, split_kms=20.0, width_kms=3.0)

Model cube and `low` / `main` / `high` layer cubes (x, y, v; Float32) from ROHSA
parameters `[amp, mu, sigma]` per Gaussian (channel units, zero-based centres),
using `reference` (channel units, per pixel) as the zero of the velocity offset.
"""
function build_layer_cubes(params, nv::Int, channel_kms::Real; reference,
                           split_kms::Real=20.0, width_kms::Real=3.0)
    n_gauss = size(params, 1) ÷ 3
    nx, ny = size(params, 2), size(params, 3)
    size(reference) == (nx, ny) || throw(ArgumentError("reference must be $nx × $ny"))
    cubes = Dict(name => zeros(Float32, nx, ny, nv) for name in (:model, ROHSA_LAYER_NAMES...))
    model, low, main, high = cubes[:model], cubes[:low], cubes[:main], cubes[:high]
    Threads.@threads for y in 1:ny
        @inbounds for x in 1:nx
            ref = reference[x, y]
            for g in 1:n_gauss
                amp = params[3g - 2, x, y]
                amp > 0 || continue
                mu = params[3g - 1, x, y]
                sigma = max(params[3g, x, y], 1e-6)
                offset = isfinite(ref) ? (mu - ref) * channel_kms : 0.0
                wl, wm, wh = velocity_layer_weights(offset; split_kms, width_kms)
                for k in 1:nv
                    v = Float32(amp * exp(-((k - 1 - mu)^2) / (2sigma^2)))
                    model[x, y, k] += v
                    low[x, y, k] += Float32(wl) * v
                    main[x, y, k] += Float32(wm) * v
                    high[x, y, k] += Float32(wh) * v
                end
            end
        end
    end
    return cubes
end
