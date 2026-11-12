# Removing the disk rotation from a spectral cube before VGT.
#
# In thin velocity channels the large-scale velocity field of a rotating disk
# carves edges into the emission (the channel picks out a band of sky where the
# disk velocity matches), and VGT reads those edges as turbulence. Shifting
# every spectrum by the local disk velocity makes the disk velocity the same
# everywhere, so channel maps no longer follow the rotation pattern.
#
# The shift is a fractional delay applied in Fourier space. It is all-pass, so
# the noise keeps the same variance and spectrum at every pixel whatever the
# fractional part of the shift; linear interpolation would smooth the noise by
# an amount that varies with the fractional shift and so imprint the contours of
# the reference velocity map onto the cube. Integer shifts are exact copies.

"""
    derotation_offset_grid(reference, nv) -> (omin, nout)

Output channel `j` holds velocity offset `omin + j - 1` channels from the
reference. The grid covers every input channel of every pixel, so no emission
falls off the ends. The reference is clamped to the band (see `clamp_reference`).
"""
function derotation_offset_grid(reference::AbstractMatrix, nv::Integer)
    finite = [clamp_reference(r, nv) for r in reference if isfinite(r)]
    isempty(finite) && throw(ArgumentError("reference has no finite values"))
    omin = floor(Int, -maximum(finite))
    omax = ceil(Int, nv - 1 - minimum(finite))
    return omin, omax - omin + 1
end

# Far from any emission the smoothed centroid is a ratio of two near-zero sums
# and can sit thousands of channels outside the band. The disk velocity of a
# pixel with emission lies inside the band, so clamping changes nothing where
# there is gas and keeps the output grid at most 2nv - 1 channels long.
clamp_reference(r::Real, nv::Integer) = clamp(Float64(r), 0.0, Float64(nv - 1))

"""
    derotate_cube(cube, reference; margin=16) -> (out, omin)

Shift every spectrum of `cube` (x, y, v) so that its velocity is measured from
`reference[x, y]` (channel units, zero-based like the ROHSA centres). Output is
Float32 on the offset grid of `derotation_offset_grid`. Spectra that are all NaN
(or have a non-finite reference) stay NaN; isolated NaN channels count as zero.

Output channels that fall outside the observed band of that pixel are returned
as NaN; `fill_padding!` then fills them (zero for noise-free model cubes,
synthetic noise for data). Where the band ends depends on the local disk
velocity, so a zero fill in a noisy cube would put a sharp noise-to-zero edge
along the disk velocity contours in the outer channels, and VGT would read that
edge as structure aligned with the rotation (exactly the imprint the
de-rotation is meant to remove). NaN cannot stay either: VGT weights each
window by the intensity of its centre pixel, so a NaN centre would poison that
pixel's Stokes sums.

`integer_shift=true` rounds every shift to a whole channel (exact copies, at
most half a channel of misalignment). Use it for cubes whose noise is not
Gaussian, such as our data with negative values clipped: a fractional shift
keeps the noise variance but changes its distribution with the fractional part
of the shift (exact clipped noise at 0, a smeared mix at 0.5), and that
fractional part repeats every channel of disk velocity, so the share of pixels
above the VGT threshold would vary in bands along the disk velocity contours.
"""
function derotate_cube(cube::AbstractArray{<:Real, 3}, reference::AbstractMatrix; margin::Int=16,
                       integer_shift::Bool=false)
    nx, ny, nv = size(cube)
    size(reference) == (nx, ny) || throw(ArgumentError("reference must be $nx × $ny"))
    omin, nout = derotation_offset_grid(reference, nv)
    len = nextprod((2, 3, 5), nout + 2margin)
    out = fill(NaN32, nx, ny, nout)
    nfreq = len ÷ 2 + 1
    FFTW.set_num_threads(1)
    nbuf = Threads.maxthreadid()
    bufs = [zeros(Float64, len) for _ in 1:nbuf]
    specs = [zeros(ComplexF64, nfreq) for _ in 1:nbuf]
    fwd = [plan_rfft(bufs[i]) for i in 1:nbuf]
    inv = [plan_irfft(specs[i], len) for i in 1:nbuf]
    Threads.@threads :static for y in 1:ny
        tid = Threads.threadid()
        buf, spec = bufs[tid], specs[tid]
        @inbounds for x in 1:nx
            isfinite(reference[x, y]) || continue
            ref = clamp_reference(reference[x, y], nv)
            integer_shift && (ref = round(ref))
            any(isfinite, view(cube, x, y, :)) || continue
            # Output j (1-based) samples input position (j - 1) + ref + omin,
            # i.e. the input delayed by d = -(ref + omin) >= 0 channels.
            d = -(ref + omin)
            di = floor(Int, d)
            df = d - di
            fill!(buf, 0.0)
            for k in 1:nv
                v = cube[x, y, k]
                buf[margin + di + k] = isfinite(v) ? Float64(v) : 0.0
            end
            if df != 0
                mul!(spec, fwd[tid], buf)
                for f in 1:nfreq
                    phase = 2pi * (f - 1) * df / len
                    if 2(f - 1) == len   # Nyquist bin of an even length must stay real
                        spec[f] *= cos(phase)
                    else
                        spec[f] *= cis(-phase)
                    end
                end
                mul!(buf, inv[tid], spec)
            end
            for j in 1:nout
                c = (j - 1) - d          # input position (zero-based) of this output channel
                out[x, y, j] = 0 <= c <= nv - 1 ? Float32(buf[margin + j]) : NaN32
            end
        end
    end
    return out, omin
end

"""
    derotated_header(header, nout, omin)

Header of a de-rotated cube: the spatial WCS is copied unchanged; the velocity
axis is the offset from the local disk velocity, in the input's units and
direction.
"""
function derotated_header(header, nout::Integer, omin::Integer)
    out = Dict{String, Any}(String(k) => header[k] for k in keys(header) if !(String(k) in ("COMMENT", "HISTORY", "")))
    cdelt = Float64(header["CDELT3"])
    out["NAXIS3"] = nout
    out["CRPIX3"] = 1.0
    out["CRVAL3"] = omin * cdelt
    out["CTYPE3"] = "VELO"   # standard linear velocity type (WCSLIB rejects invented ones)
    out["DEROTATE"] = "velocity = offset from the local disk velocity"
    return out
end

"""
    noise_sigma_map(cube; smooth_px=15.0) -> Matrix{Float64}

Per-pixel noise of a data cube: robust estimate from the negative samples of
each spectrum, smoothed spatially (the per-spectrum estimate from ~100 channels
scatters by ~10%). Pixels without data get the median.
"""
function noise_sigma_map(cube::AbstractArray{<:Real, 3}; smooth_px::Real=15.0)
    nx, ny, _ = size(cube)
    sigma = fill(NaN, nx, ny)
    Threads.@threads for y in 1:ny
        for x in 1:nx
            sigma[x, y] = robust_negative_noise(view(cube, x, y, :); min_count=20)
        end
    end
    finite = filter(isfinite, sigma)
    isempty(finite) && throw(ArgumentError("no negative samples to estimate the noise from (clipped cube?); pass the noise level instead"))
    fill_value = median(finite)
    sigma = map(s -> isfinite(s) ? s : fill_value, sigma)
    return imfilter(sigma, Kernel.gaussian((Float64(smooth_px), Float64(smooth_px))), "reflect")
end

"""
    fill_padding!(out, observed; sigma=nothing, beam_sigma_px=1.27, clip=false, seed=1)

Replace the NaN padding of `derotate_cube` output. `observed[x, y]` marks pixels
that have data (others stay NaN). With `sigma === nothing` the padding becomes
zero (noise-free model cubes). With a noise map it becomes Gaussian noise of
that per-pixel level, smoothed spatially with the beam (`beam_sigma_px`, 1′ FWHM
on 20″ pixels = 1.27 px) and renormalised, so the padded part of each channel
has the same amplitude and spatial correlation as the observed noise and there
is no edge for VGT to find. `clip=true` sets negative noise to zero, for cubes
whose negative values were clipped in preprocessing (the padding must look like
the observed noise, including that clipping).
"""
function fill_padding!(out::AbstractArray{Float32, 3}, observed::AbstractMatrix{Bool};
                       sigma=nothing, beam_sigma_px::Real=1.27, clip::Bool=false, seed::Integer=1)
    nx, ny, nout = size(out)
    if sigma === nothing
        @inbounds for k in 1:nout, y in 1:ny, x in 1:nx
            if observed[x, y] && isnan(out[x, y, k])
                out[x, y, k] = 0f0
            end
        end
        return out
    end
    kernel = Kernel.gaussian((Float64(beam_sigma_px), Float64(beam_sigma_px)))
    Threads.@threads for k in 1:nout
        plane = view(out, :, :, k)
        any(i -> observed[i] && isnan(plane[i]), eachindex(plane)) || continue
        rng = Random.Xoshiro(seed + k)
        white = imfilter(randn(rng, nx, ny), kernel, "circular")
        white ./= std(white)
        @inbounds for i in eachindex(plane)
            if observed[i] && isnan(plane[i])
                v = sigma[i] * white[i]
                plane[i] = Float32(clip ? max(v, 0.0) : v)
            end
        end
    end
    return out
end

"""
    sanitize_reference(reference, column, nv; lo_frac=0.0005, hi_frac=0.002, smooth_px=5.0)

Reference velocity that can be trusted for shifting. The smoothed model centroid
is meaningless where the model has (almost) no emission: there it is a ratio of
two near-zero sums and jumps by hundreds of channels between neighbours.
Shifting those spectra by such amounts would create structure along the jumps.
Where the model column `column` (sum over channels) is below `lo_frac` of its
99th percentile the reference becomes the median disk velocity of the trusted
area; above `hi_frac` it is kept; in between the two are blended, so the shift
never jumps. On the LMC 2048² window the reference stays within 30–115
channels wherever the column exceeds 0.1% of its 99th percentile (83% of the
pixels) and reaches −4224…+7936 channels only where there is no emission.
"""
function sanitize_reference(reference::AbstractMatrix, column::AbstractMatrix, nv::Integer;
                            lo_frac::Real=0.0005, hi_frac::Real=0.002, smooth_px::Real=5.0)
    size(reference) == size(column) || throw(ArgumentError("reference and column must have the same size"))
    top = quantile(filter(isfinite, vec(column)), 0.99)
    weight = map(c -> isfinite(c) ? clamp((c / top - lo_frac) / (hi_frac - lo_frac), 0.0, 1.0) : 0.0, column)
    weight = clamp.(imfilter(weight, Kernel.gaussian((Float64(smooth_px), Float64(smooth_px))), "replicate"), 0.0, 1.0)
    trusted = [weight[i] >= 0.999 && isfinite(reference[i]) for i in eachindex(reference)]
    any(trusted) || throw(ArgumentError("no pixel with enough emission to trust the reference"))
    r0 = median(clamp_reference.(reference[trusted], nv))
    out = similar(weight)
    for i in eachindex(out)
        r = isfinite(reference[i]) ? clamp_reference(reference[i], nv) : r0
        out[i] = weight[i] * r + (1 - weight[i]) * r0
    end
    return out, weight
end
