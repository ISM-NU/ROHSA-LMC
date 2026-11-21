using DrWatson
Base.exit_on_sigint(false)
quickactivate(@__DIR__, "RohsaLMC")

using FITSIO
using Printf
using Statistics
using VGT
# Reliable channel range per pixel for ROHSA's `channel_range_file` (with missing_data = "skip").
#
# Why: in a mosaic of observing blocks with different channel ranges, the end-of-band channels that only one of two
# overlapping blocks covers can be averaged with the empty block: signal and noise are scaled by a weight w that falls
# towards the seam. Those voxels are not 0, so the live range keeps them. This tool measures w per channel and tile
# from the noise (rohsa/inputs.jl attenuation_weight) and keeps, per pixel, the outermost live channels with
# w >= --wmin (reliable_channel_ranges). No 1/w renormalisation: attenuated channels are dropped, not rescaled.
#
#   julia --project=. --threads=32 rohsa/cli/channel_ranges.jl --cube <cube.fits> --out <ranges.fits> \
#       [--x0 1 --y0 1 --size 0] [--tile 16] [--smooth 5] [--wmin 0.8] [--ref K]
#
# --x0/--y0/--size: crop as in the ROHSA configs (crop_x0, crop_y0, crop_size; 1-based; size 0 = whole image).
# Output FITS (Julia x-y order, 1-based channels; first > last = no data): FIRST, LAST (reliable range), LIVEFIRST,
# LIVELAST (first/last channel neither 0 nor NaN), W (tile x tile x channel weights; header REF = reference noise K,
# TILE). Memory: the crop as Float32 (2.7 GB for 2048^2 x 159) plus about 1 GB.

Base.include(@__MODULE__, joinpath(@__DIR__, "..", "inputs.jl"))

function main(args)
    o = Dict("--cube" => "", "--out" => "", "--x0" => "1", "--y0" => "1", "--size" => "0", "--tile" => "16",
             "--smooth" => "5", "--wmin" => "0.8", "--ref" => "")
    i = 1
    while i <= length(args)
        haskey(o, args[i]) || throw(ArgumentError("unknown option $(args[i])"))
        o[args[i]] = args[i + 1]
        i += 2
    end
    (isempty(o["--cube"]) || isempty(o["--out"])) && throw(ArgumentError("--cube and --out are required"))
    x0, y0, n = parse(Int, o["--x0"]), parse(Int, o["--y0"]), parse(Int, o["--size"])
    tile = parse(Int, o["--tile"])
    f = FITS(o["--cube"])
    cube = try
        nx, ny, _ = size(f[1])
        xs = n > 0 ? (x0:x0 + n - 1) : (1:nx)
        ys = n > 0 ? (y0:y0 + n - 1) : (1:ny)
        Float32.(read(f[1], xs, ys, :))
    finally
        close(f)
    end
    first, last = live_channel_ranges(cube)
    ref = isempty(o["--ref"]) ? nothing : parse(Float64, o["--ref"])
    t = @elapsed w, r = attenuation_weight(cube, first, last; tile=tile,
                                           smooth=parse(Int, o["--smooth"]), ref=ref)
    wmin = parse(Float64, o["--wmin"])
    rf, rl = reliable_channel_ranges(first, last, w; tile=tile, wmin=wmin)
    live = first .<= last
    cut = (rf .- first) .+ (last .- rl)
    @printf("reference noise %.4f K (%s), weights in %.0f s\n", r, ref === nothing ? "mode over tiles and channels" : "given", t)
    @printf("pixels with data %d; with channels cut %d (%.3f); median channels cut where any %.0f; no reliable channel %d\n",
            count(live), count(live .& (cut .> 0)), count(live .& (cut .> 0)) / count(live),
            count(live .& (cut .> 0)) > 0 ? median(cut[live .& (cut .> 0)]) : 0, count(live .& (rf .> rl)))
    FITS(o["--out"], "w") do g
        hdr = FITSHeader(["CUBE", "X0", "Y0", "SIZE", "WMIN"], [basename(o["--cube"]), x0, y0, n, wmin],
                         ["", "crop x0 (1-based)", "crop y0", "crop size (0 = whole)", "minimum weight kept"])
        write(g, Int16.(rf); header=hdr, name="FIRST")
        write(g, Int16.(rl); name="LAST")
        write(g, Int16.(first); name="LIVEFIRST")
        write(g, Int16.(last); name="LIVELAST")
        whdr = FITSHeader(["REF", "TILE"], [r, tile], ["reference noise (K)", "tile size (px)"])
        write(g, Float32.(w); header=whdr, name="W")
    end
    println("wrote ", o["--out"])
end

main(ARGS)
