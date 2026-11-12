using DrWatson
Base.exit_on_sigint(false)
quickactivate(@__DIR__, "RohsaLMC")

using Statistics
using TOML
using Astronomy
using RohsaLMC
# De-rotate cubes by the local disk velocity (src/science/derotate.jl) and write
# a VGT config for each, copied from an existing config with only the input and
# output paths changed, so the de-rotated run differs from the original in
# nothing but the de-rotation.
#
#   julia --project=. rohsa/cli/derotate_cube.jl <reference.fits> <out_dir> <config_dir> \
#       [--noise-pad data] [--beam-fwhm-px 3] [--column model_moment0.fits] [--integer-shift data] name=cube.fits ...
#
# <reference.fits>: disk velocity in channel units (reference_velocity_channels.fits
# from rohsa/cli/layer_cubes.jl). <config_dir> holds vgt_<name>.toml for every name.
# Channels outside a pixel's observed band are filled with beam-correlated noise
# for the cubes named in --noise-pad (observed data) and with zero otherwise
# (noise-free ROHSA model and layer cubes); see fill_padding!. For a cube with
# clipped negatives the noise level is the [noise] noise of its VGT config and
# the padding noise is clipped too. --column (model integrated intensity) makes
# the shift constant where the model has no emission (sanitize_reference).
# --integer-shift lists cubes shifted by whole channels (clipped, non-Gaussian
# noise; see derotate_cube).

const USAGE = "Usage: julia --project=. rohsa/cli/derotate_cube.jl <reference.fits> <out_dir> <config_dir> [--noise-pad data] [--beam-fwhm-px 3] [--column model_moment0.fits] [--integer-shift data] name=cube.fits ..."

function main(args=ARGS)
    length(args) >= 4 || throw(ArgumentError(USAGE))
    reference_path, out_dir, config_dir = abspath.(args[1:3])
    rest = collect(args[4:end])
    noise_pad = ["data"]
    beam_fwhm_px = 3.0
    column_path = ""
    integer_names = String[]
    while !isempty(rest) && startswith(rest[1], "--")
        length(rest) >= 2 || throw(ArgumentError(USAGE))
        opt, val = popfirst!(rest), popfirst!(rest)
        if opt == "--noise-pad"
            noise_pad = split(val, r"[,:]")   # ":" also, because sbatch --export splits on ","
        elseif opt == "--beam-fwhm-px"
            beam_fwhm_px = parse(Float64, val)
        elseif opt == "--column"
            column_path = abspath(val)
        elseif opt == "--integer-shift"
            integer_names = split(val, r"[,:]")
        else
            throw(ArgumentError("Unknown option $opt. $USAGE"))
        end
    end
    inputs = map(rest) do arg
        occursin('=', arg) || throw(ArgumentError(USAGE))
        name, path = split(arg, '='; limit=2)
        String(name) => abspath(path)
    end
    reference = Float64.(read_fits_image(reference_path))
    mkpath(out_dir)
    trace = Dict{String, Any}("reference" => reference_path, "cubes" => Dict{String, Any}())
    if !isempty(column_path)
        nv_first = Int(read_fits_header(last(first(inputs)))["NAXIS3"])
        reference, weight = RohsaLMC.sanitize_reference(reference, Float64.(read_fits_image(column_path)), nv_first)
        write_fits(joinpath(out_dir, "reference_velocity_channels_sanitized.fits"), Float32.(reference); header=nothing)
        write_fits(joinpath(out_dir, "reference_trust_weight.fits"), Float32.(weight); header=nothing)
        trace["column"] = column_path
        trace["trusted_fraction"] = count(>=(0.999), weight) / length(weight)
        println("reference sanitized with $column_path: trusted fraction $(trace["trusted_fraction"])")
    end
    for (name, path) in inputs
        template = joinpath(config_dir, "vgt_$(name).toml")
        isfile(template) || throw(ArgumentError("Missing VGT config $template"))
        header = read_fits_header(path)
        cube = read_fits_image(path)
        t = @elapsed out, omin = RohsaLMC.derotate_cube(cube, reference; integer_shift=name in integer_names)
        flux_in = sum(v -> isfinite(v) ? Float64(v) : 0.0, cube)
        flux_out = sum(v -> isfinite(v) ? Float64(v) : 0.0, out)
        observed = dropdims(any(isfinite, cube; dims=3); dims=3)
        if name in noise_pad
            # Our preprocessing clips negative values; then the noise level cannot be
            # measured from the negative tail, so take the one VGT itself uses.
            clipped = minimum(v -> isfinite(v) ? v : Inf, cube) >= 0
            sigma = clipped ? fill(Float64(TOML.parsefile(template)["noise"]["noise"]), size(observed)) :
                              RohsaLMC.noise_sigma_map(cube)
            println("$(name): padding with beam-correlated noise, clipped=$(clipped), sigma median $(median(sigma)) range $(extrema(sigma))")
            RohsaLMC.fill_padding!(out, observed; sigma, clip=clipped, beam_sigma_px=beam_fwhm_px / sqrt(8log(2)))
        else
            RohsaLMC.fill_padding!(out, observed)
        end
        cube = nothing
        stem = splitext(basename(path))[1]
        out_path = joinpath(out_dir, "$(stem)_derot.fits")
        write_fits(out_path, out; header=RohsaLMC.derotated_header(header, size(out, 3), omin))
        println("$(name): $(size(out)) offsets $(omin)..$(omin + size(out, 3) - 1) channels, ",
                "flux ratio $(flux_out / flux_in), $(round(t; digits=1)) s -> $out_path")
        cfg = TOML.parsefile(template)
        cfg["input"]["input_file"] = out_path
        cfg["input"]["output_dir"] = joinpath(out_dir, "vgt", name)
        open(io -> TOML.print(io, cfg; sorted=true), joinpath(out_dir, "vgt_$(name).toml"), "w")
        trace["cubes"][name] = Dict("input" => path, "output" => out_path, "omin" => omin,
                                    "nout" => size(out, 3), "flux_ratio" => flux_out / flux_in,
                                    "padding" => name in noise_pad ? "beam-correlated noise (clipped if the input is clipped)" : "zero",
                                    "shift" => name in integer_names ? "integer channels" : "fractional (Fourier)")
        out = nothing
        GC.gc()
    end
    open(io -> TOML.print(io, trace; sorted=true), joinpath(out_dir, "derotate_trace.toml"), "w")
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end
