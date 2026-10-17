using FITSIO
using TOML
using Astronomy

include(joinpath(@__DIR__, "..", "..", "memory.jl"))
include(joinpath(@__DIR__, "..", "..", "configs.jl"))

const IRGETAS_ZEN4_USABLE_GIB = 385_974 / 1024
const IRGETAS_H100_USABLE_GIB = 770_000 / 1024
const IRGETAS_HEADROOM = 0.90

function irgetas_partition(peak_gib::Real)
    peak_gib <= IRGETAS_ZEN4_USABLE_GIB * IRGETAS_HEADROOM && return :ZEN4
    peak_gib <= IRGETAS_H100_USABLE_GIB * IRGETAS_HEADROOM && return :H100
    return :none
end

function inspect_config(config_path::AbstractString)
    opts = load_rohsa_config(config_path)
    f = FITS(opts["input"], "r")
    nx, ny, nv = try
        size(f[1])
    finally
        close(f)
    end
    crop = parse(Int, opts["crop-size"])
    stride = parse(Int, opts["spatial-stride"])
    source_nx = crop > 0 ? crop : nx
    source_ny = crop > 0 ? crop : ny
    out_nx = cld(source_nx, stride)
    out_ny = cld(source_ny, stride)
    ns = parse.(Int, split(opts["n-gauss"], ','))
    estimate = rohsa_memory_estimate(out_nx, out_ny, nv, maximum(ns), parse(Int, opts["m"]))
    return (; opts, nx=out_nx, ny=out_ny, nv, estimate,
            partition=irgetas_partition(estimate.peak_gib))
end

function main(args)
    length(args) == 1 || error("Usage: julia --project=. rohsa/hpc/irgetas/preflight.jl CONFIG.toml")
    result = inspect_config(args[1])
    e = result.estimate
    println("ROHSA output grid: $(result.nx) x $(result.ny) x $(result.nv)")
    println("ROHSA padded grid: $(e.padded_side) x $(e.padded_side)")
    println("Conservative peak estimate: $(round(e.peak_gib; digits=1)) GiB")
    if result.partition == :ZEN4
        println("Irgetas result: fits the conservative ZEN4 safety limit.")
    elseif result.partition == :H100
        println("Irgetas result: H100 memory is required; ZEN4 is not safe.")
        println("ROHSA does not use the GPUs. Ask the NU HPC team before a CPU-only H100 submission.")
    else
        println(stderr, "Irgetas result: configuration exceeds the conservative H100 safety limit.")
        return 2
    end
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main(ARGS))
end
