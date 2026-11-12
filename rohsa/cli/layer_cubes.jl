using DrWatson
Base.exit_on_sigint(false)
quickactivate(@__DIR__, "RohsaLMC")

using TOML
using VGT
using RohsaLMC
# Build velocity-layer cubes (low / main / high) and the full model cube from one
# ROHSA region directory, plus a VGT config per cube and one for the data crop
# itself (the control: VGT on data vs on the model shows what the
# decomposition alone changes).
#
#   julia --project=. rohsa/cli/layer_cubes.jl <region_dir> <n_gauss> <out_dir> \
#       [--split-kms 20] [--width-kms 3] [--ref-sigma-px 30] [--noise K] [--ny-chunks N]
#       [--split-by velocity|width] [--cold-min-amp K]
#
# --split-by width writes cold / unstable / warm cubes instead (Gaussians by
# Doppler temperature: < 500 K, 500-5000 K, >= 5000 K); --cold-min-amp drops
# narrow Gaussians fainter than K (fitted noise). --cold-max-fwhm / --warm-min-fwhm
# move those edges to line widths in km/s; giving both the same value makes a
# two-way split (narrow / broad) with an empty unstable cube.

const USAGE = "Usage: julia --project=. rohsa/cli/layer_cubes.jl <region_dir> <n_gauss> <out_dir> [--split-kms 20] [--width-kms 3] [--ref-sigma-px 30] [--noise K] [--ny-chunks N] [--split-by velocity|width] [--cold-min-amp K] [--phase-log-width W] [--cold-max-fwhm KMS] [--warm-min-fwhm KMS]"

function parse_layer_args(args)
    length(args) >= 3 || throw(ArgumentError(USAGE))
    opts = Dict("--split-kms" => "20", "--width-kms" => "3", "--ref-sigma-px" => "30",
                "--noise" => "", "--ny-chunks" => "16",
                "--split-by" => "velocity", "--cold-min-amp" => "0", "--phase-log-width" => "0.2",
                "--cold-max-fwhm" => "", "--warm-min-fwhm" => "")
    rest = args[4:end]
    iseven(length(rest)) || throw(ArgumentError(USAGE))
    for i in 1:2:length(rest)
        haskey(opts, rest[i]) || throw(ArgumentError("Unknown option $(rest[i]). $USAGE"))
        opts[rest[i]] = rest[i + 1]
    end
    return (; region_dir=abspath(args[1]), n_gauss=parse(Int, args[2]), out_dir=abspath(args[3]),
            split_kms=parse(Float64, opts["--split-kms"]), width_kms=parse(Float64, opts["--width-kms"]),
            ref_sigma_px=parse(Float64, opts["--ref-sigma-px"]),
            noise=isempty(opts["--noise"]) ? nothing : parse(Float64, opts["--noise"]),
            ny_chunks=parse(Int, opts["--ny-chunks"]), split_by=opts["--split-by"],
            cold_min_amp=parse(Float64, opts["--cold-min-amp"]), phase_log_width=parse(Float64, opts["--phase-log-width"]),
            # Phase edges in K (defaults 500 / 5000 K), or from line widths when given in km/s.
            cold_max_k=isempty(opts["--cold-max-fwhm"]) ? 500.0 :
                RohsaLMC.fwhm_doppler_temperature(parse(Float64, opts["--cold-max-fwhm"])),
            warm_min_k=isempty(opts["--warm-min-fwhm"]) ? 5000.0 :
                RohsaLMC.fwhm_doppler_temperature(parse(Float64, opts["--warm-min-fwhm"])))
end

# The ROHSA wrapper writes one scalar noise value to every pixel of rms_map.dat.
function rohsa_scalar_noise(region_dir)
    open(joinpath(region_dir, "rms_map.dat")) do io
        readline(io)
        return parse(Float64, split(readline(io))[3])
    end
end

function write_vgt_config(path, input_file, output_dir, noise, ny_chunks)
    cfg = TOML.parsefile(joinpath(dirname(@__DIR__), "..", "config.example.toml"))
    cfg["input"]["input_file"] = input_file
    cfg["input"]["output_dir"] = output_dir
    cfg["input"]["keep_intermediate_chunks"] = false
    cfg["noise"]["noise"] = noise
    cfg["execution"]["ny_chunks"] = ny_chunks
    open(io -> TOML.print(io, cfg; sorted=true), path, "w")
    return path
end

function main(args=ARGS)
    opt = parse_layer_args(args)
    region = basename(rstrip(opt.region_dir, '/'))
    data_path = joinpath(opt.region_dir, "$(region).fits")
    rohsa_out = joinpath(opt.region_dir, "rohsa_n$(opt.n_gauss).dat")
    isfile(data_path) || throw(ArgumentError("Missing data crop $data_path"))
    isfile(rohsa_out) || throw(ArgumentError("Missing ROHSA output $rohsa_out"))
    mkpath(opt.out_dir)

    header = read_fits_header(data_path)
    nx, ny, nv = (Int(header["NAXIS$i"]) for i in 1:3)
    channel_kms = abs(Float64(header["CDELT3"])) / 1000
    noise = something(opt.noise, rohsa_scalar_noise(opt.region_dir))
    println("region=$region grid=$nx×$ny×$nv channel=$(channel_kms) km/s noise=$noise K")

    params = RohsaLMC.read_rohsa_parameters(rohsa_out, opt.n_gauss, nx, ny)
    opt.split_by in ("velocity", "width") || throw(ArgumentError("--split-by must be velocity or width"))
    reference = RohsaLMC.reference_centroid_channels(params; sigma_pixels=opt.ref_sigma_px)
    cubes = opt.split_by == "width" ?
        RohsaLMC.build_phase_cubes(params, nv, channel_kms; cold_max_k=opt.cold_max_k, warm_min_k=opt.warm_min_k,
                                    cold_min_amp=opt.cold_min_amp, log_width=opt.phase_log_width) :
        RohsaLMC.build_layer_cubes(params, nv, channel_kms; reference,
                                    split_kms=opt.split_kms, width_kms=opt.width_kms)

    model_total = sum(Float64, cubes[:model])
    configs = String[]
    write_fits(joinpath(opt.out_dir, "reference_velocity_channels.fits"), Float32.(reference); header=nothing)
    for (name, cube) in cubes
        path = joinpath(opt.out_dir, "$(region)_$(name).fits")
        write_fits(path, cube; header)
        share = sum(Float64, cube) / model_total
        println("$(name): $(round(100share; digits=1))% of model emission -> $path")
        if name != :model
            fraction = dropdims(sum(cube; dims=3); dims=3) ./ dropdims(sum(cubes[:model]; dims=3); dims=3)
            write_fits(joinpath(opt.out_dir, "$(region)_$(name)_fraction.fits"), Float32.(fraction); header=nothing)
        end
        push!(configs, write_vgt_config(joinpath(opt.out_dir, "vgt_$(name).toml"), path,
                                        joinpath(opt.out_dir, "vgt", String(name)), noise, opt.ny_chunks))
    end
    push!(configs, write_vgt_config(joinpath(opt.out_dir, "vgt_data.toml"), data_path,
                                    joinpath(opt.out_dir, "vgt", "data"), noise, opt.ny_chunks))
    open(joinpath(opt.out_dir, "layers_trace.toml"), "w") do io
        TOML.print(io, Dict("region_dir" => opt.region_dir, "n_gauss" => opt.n_gauss,
                            "split_kms" => opt.split_kms, "width_kms" => opt.width_kms,
                            "ref_sigma_px" => opt.ref_sigma_px, "noise" => noise,
                            "split_by" => opt.split_by, "cold_min_amp" => opt.cold_min_amp, "phase_log_width" => opt.phase_log_width,
                            "cold_max_k" => opt.cold_max_k, "warm_min_k" => opt.warm_min_k,
                            "vgt_configs" => configs); sorted=true)
    end
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end
