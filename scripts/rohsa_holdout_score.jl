using DrWatson
quickactivate(@__DIR__, "RohsaLMC")

using FITSIO
using Statistics
using RohsaLMC: read_rohsa_parameters

# Score a ROHSA run by prediction: fit with a held-out mask (rohsa config holdout_mask, noise 0 there, so the
# full-resolution fit sees no data in those pixels), then compare the model with the data in the held-out block
# centres. Extra freedom (more Gaussians, more iterations, weaker smoothing) always lowers the in-sample chi2,
# also when it only fits noise; it should not improve the prediction of spectra the fit never saw.
#
# On a null cube the noise-free truth is known; --truth-* gives it, and the script also reports the model's error
# against the truth, in-sample and held out. That is the check that held-out scores rank configurations the way
# the true error does (if they do not, they cannot be used to choose configurations on the data).
#
#   julia --project=. scripts/rohsa_holdout_score.jl --rohsa run/crop/rohsa_n16.dat --n-gauss 16 \
#       --cube run/crop/crop.fits --mask holdout.fits --sigma 0.236 --out score.csv [--core-dist 7] \
#       [--truth-dat region/rohsa_n16.dat --truth-n 16 --truth-size 2048 --truth-offset 1536,1024 \
#        --truth-broaden-kms 3 --dv 1]
#
# --core-dist d: held-out pixels at Chebyshev distance >= d from any fitted pixel are scored (7 in a 16-px block
# = its central 4x4), so the smoothing cannot simply copy a fitted neighbour. Voxels that are exactly 0 in the
# cube (its blank channels) are skipped. --truth-broaden-kms: the null v2 rule (make_rohsa_null_cube.jl
# --broaden): truth components narrower than this are widened to it with a*sigma kept. CNM share = flux in
# Gaussians with sigma < 2.03 km/s; on a null it is invented.

const CNM_KMS = 2.03

function parse_args(args)
    o = Dict{String,String}("core-dist" => "7", "dv" => "1", "truth-dat" => "", "truth-broaden-kms" => "0")
    i = 1
    while i <= length(args)
        startswith(args[i], "--") && i < length(args) || throw(ArgumentError("Bad argument: $(args[i])"))
        o[args[i][3:end]] = args[i + 1]
        i += 2
    end
    for k in ("rohsa", "n-gauss", "cube", "mask", "sigma", "out")
        haskey(o, k) || throw(ArgumentError("--$k is required"))
    end
    return o
end

gauss_spectrum!(m, p, x, y, n) = begin
    fill!(m, 0.0)
    for g in 1:n
        a, mu, s = p[3g - 2, x, y], p[3g - 1, x, y], p[3g, x, y]
        (a == 0 || !(s > 0)) && continue
        for k in eachindex(m)
            m[k] += a * exp(-0.5 * ((k - 1 - mu) / s)^2)
        end
    end
    m
end

function broaden!(p, n, min_ch)
    min_ch > 0 || return p
    for y in axes(p, 3), x in axes(p, 2), g in 1:n
        a, s = p[3g - 2, x, y], p[3g, x, y]
        if s > 0 && s < min_ch
            p[3g - 2, x, y] = a * s / min_ch
            p[3g, x, y] = min_ch
        end
    end
    p
end

# Chebyshev distance of each held-out pixel to the nearest fitted pixel, capped at cap.
function holdout_distance(held, cap)
    nx, ny = size(held)
    d = zeros(Int, nx, ny)
    for y in 1:ny, x in 1:nx
        held[x, y] || continue
        r = 1
        while r < cap
            found = false
            for yy in max(1, y - r):min(ny, y + r), xx in max(1, x - r):min(nx, x + r)
                max(abs(xx - x), abs(yy - y)) == r && !held[xx, yy] && (found = true; break)
            end
            found && break
            r += 1
        end
        d[x, y] = r
    end
    d
end

function main(args)
    o = parse_args(args)
    n = parse(Int, o["n-gauss"]); sigma = parse(Float64, o["sigma"]); dv = parse(Float64, o["dv"])
    cube = Float64.(FITS(f -> read(f[1]), o["cube"]))  # (x, y, v)
    nx, ny, nv = size(cube)
    held = FITS(f -> read(f[1]), o["mask"]) .!= 0
    size(held) == (nx, ny) || throw(ArgumentError("mask is $(size(held)), cube is $((nx, ny))"))
    p = read_rohsa_parameters(o["rohsa"], n, nx, ny)
    dist = holdout_distance(held, parse(Int, o["core-dist"]))
    core = dist .>= parse(Int, o["core-dist"])
    truth = nothing
    if !isempty(o["truth-dat"])
        tn = parse(Int, o["truth-n"]); ts = parse(Int, o["truth-size"])
        ox, oy = parse.(Int, split(o["truth-offset"], ","))
        tp = read_rohsa_parameters(o["truth-dat"], tn, ts, ts)[:, ox + 1:ox + nx, oy + 1:oy + ny]
        truth = (p=broaden!(tp, tn, parse(Float64, o["truth-broaden-kms"]) / dv), n=tn)
    end
    sets = Dict("held_core" => core, "fitted" => .!held)
    rows = Dict{String,Any}()
    m = zeros(nv); t = zeros(nv)
    for (name, sel) in sets
        r2 = 0.0; e2 = 0.0; nvox = 0; cnm = 0.0; tot = 0.0; npix = 0
        for y in 1:ny, x in 1:nx
            sel[x, y] || continue
            npix += 1
            gauss_spectrum!(m, p, x, y, n)
            truth === nothing || gauss_spectrum!(t, truth.p, x, y, truth.n)
            for k in 1:nv
                v = cube[x, y, k]
                (v == 0 || !isfinite(v)) && continue
                r2 += (v - m[k])^2; nvox += 1
                truth === nothing || (e2 += (m[k] - t[k])^2)
            end
            for g in 1:n
                a, s = p[3g - 2, x, y], p[3g, x, y]
                (a > 0 && s > 0) || continue
                tot += a * s; s * dv < CNM_KMS && (cnm += a * s)
            end
        end
        rows["$(name)_pixels"] = npix
        rows["$(name)_chi2_per_voxel"] = r2 / max(nvox, 1) / sigma^2
        rows["$(name)_cnm_share"] = cnm / max(tot, eps())
        truth === nothing || (rows["$(name)_truth_rmse_K"] = sqrt(e2 / max(nvox, 1)))
    end
    ks = sort(collect(keys(rows)))
    open(o["out"], "w") do io
        println(io, join(ks, ","))
        println(io, join([string(rows[k]) for k in ks], ","))
    end
    for k in ks
        println(k, " = ", rows[k])
    end
end

main(ARGS)
