# ROHSA config loading and parameter-file serialization.
#
# Provenance:
# - `configs/rohsa/full_cube.toml` starts from rohsa/upstream/ROHSA/ROHSApy.ipynb,
#   parameter cell around the notebook section "Generate files for ROHSA with
#   your data".

const ROHSA_CONFIG_KEYS = Dict(
    "input" => "input",
    "output_dir" => "output-dir",
    "crop_size" => "crop-size",
    "crop_x0" => "crop-x0",
    "crop_y0" => "crop-y0",
    "spatial_stride" => "spatial-stride",
    "reuse" => "reuse",
    "write_visuals" => "write-visuals",
    "derive_products" => "derive-products",
    "write_region_fits" => "write-region-fits",
    "write_cubes" => "write-cubes",
    "keep_rohsa_input" => "keep-rohsa-input",
    "preflight_only" => "preflight-only",
    "force_full_cube" => "force-full-cube",
    "memory_gb" => "memory-gb",
    "safe_memory_fraction" => "safe-memory-fraction",
    "fits_extension" => "fits-extension",
    "spectral_kernel" => "spectral-kernel",
    "clip_negative" => "clip-negative",
    "noise_map" => "noise-map",
    "holdout_mask" => "holdout-mask",
    "noise_estimate" => "noise-estimate",
    "noise_per_pixel" => "noise-per-pixel",
    "noise_channels" => "noise-channels",
    "noise_floor_quantile" => "noise-floor-quantile",
    "missing_data" => "missing-data",
    "channel_range_file" => "channel-range-file",
    "n_gauss" => "n-gauss",
    "lambda_amp" => "lambda-amp",
    "lambda_mu" => "lambda-mu",
    "lambda_sig" => "lambda-sig",
    "lambda_var_amp" => "lambda-var-amp",
    "lambda_var_mu" => "lambda-var-mu",
    "lambda_var_sig" => "lambda-var-sig",
    "amp_fact_init" => "amp-fact-init",
    "sig_init" => "sig-init",
    "lb_sig_init" => "lb-sig-init",
    "ub_sig_init" => "ub-sig-init",
    "lb_sig" => "lb-sig",
    "ub_sig" => "ub-sig",
    "regul_grad_adjoint" => "regul-grad-adjoint",
    "init_option" => "init-option",
    "maxiter_init" => "maxiter-init",
    "maxiter" => "maxiter",
    "maxiter_last" => "maxiter-last",
    "checkpoint_at" => "checkpoint-at",
    "m" => "m",
    "noise" => "rohsa-noise",
    "regul" => "regul",
    "descent" => "descent",
    "lstd" => "lstd",
    "ustd" => "ustd",
    "iprint" => "iprint",
    "iprint_init" => "iprint-init",
    "save_grid" => "save-grid",
)

const ROHSA_AUTO_CONFIG_KEYS = Dict(
    "enabled" => "auto-config-enabled",
    "profile" => "auto-config-profile",
    "max_candidates" => "auto-config-max-candidates",
    "n_gauss" => "auto-config-n-gauss",
    "lambda_amp" => "auto-config-lambda-amp",
    "lambda_mu" => "auto-config-lambda-mu",
    "lambda_sig" => "auto-config-lambda-sig",
    "lambda_var_amp" => "auto-config-lambda-var-amp",
    "lambda_var_mu" => "auto-config-lambda-var-mu",
    "lambda_var_sig" => "auto-config-lambda-var-sig",
    "small_component_fraction" => "auto-config-small-component-fraction",
    "broad_fraction_limit" => "auto-config-broad-fraction-limit",
    "sigma_bound_tolerance" => "auto-config-sigma-bound-tolerance",
    "duplicate_velocity_kms" => "auto-config-duplicate-velocity-kms",
    "duplicate_sigma_kms" => "auto-config-duplicate-sigma-kms",
    "physical_weight" => "auto-config-physical-weight",
    "stability_weight" => "auto-config-stability-weight",
)

const OPTIONAL_ROHSA_OPTIONS = Dict(
    "spatial-stride" => "1",
    "derive-products" => "true",
    "write-region-fits" => "true",
    # 3D model/component/group cubes are rebuildable from the parameter file
    # (2.6 GB each at 2048^2); the 2D maps are always written.
    "write-cubes" => "false",
    # The text copy of the input cube (~20 GB at 2048^2) is deleted after the
    # region's ROHSA runs succeed unless this is true.
    "keep-rohsa-input" => "false",
    "preflight-only" => "false",
    "force-full-cube" => "false",
    "memory-gb" => "auto",
    "safe-memory-fraction" => "0.70",
    "fits-extension" => ".fits",
    # Spectral smoothing applied to the raw crop before clipping (rohsa/inputs.jl
    # spectral_kernel_weights); fitted widths then include the kernel in quadrature.
    "spectral-kernel" => "none",
    # Zero negative voxels before ROHSA (historical default). Upstream ROHSA does
    # not require it, and on noise-dominated cubes (CO, line-free H I channels)
    # it turns zero-mean noise into a positive pedestal that ROHSA fits with
    # broad components. NaNs are always written as 0.
    "clip-negative" => "true",
    # Exact gradient of ROHSA's smoothing term at map edges (mod_array.f90 regul_grad_adjoint). The
    # original code is wrong in the two outer rows and columns of every level. On by default since
    # 2026-11-25 (same fix as antoinemarchal/ROHSA-FITS#4); set "false" to reproduce runs made before.
    # Emitted only when true, so "false" also works with ROHSA binaries that lack the switch.
    "regul-grad-adjoint" => "true",
    # Budget study (mod_array.f90): evaluation limit of the full-resolution level only (0 = maxiter), and
    # evaluation counts at which that level's parameters are written to <run>/checkpoints/<run label>_e<k>.ckpt
    # (read with read_rohsa_checkpoint). Emitted only when set, so older binaries still read the namelist.
    "maxiter-last" => "0",
    "checkpoint-at" => "",
    # Per-pixel noise for ROHSA (2D FITS, crop-sized, Julia x-y order). Empty = one level from the cube
    # (robust_channel_difference_noise). Pixels with NaN or <= 0 are written as 0, which ROHSA skips in the data
    # term at the full-resolution level (their parameters follow the smoothing term): used for held-out blocks.
    "noise-map" => "",
    # Held-out pixels (2D FITS, crop-sized, nonzero = held out): their noise is set to 0, so ROHSA fits them from
    # the smoothing term only (full-resolution level; coarse levels still average them in). Scoring the fit on
    # them measures how well a configuration predicts spectra it did not fit.
    "holdout-mask" => "",
    # How the noise level for rms_map.dat is measured. "channel_difference" (default, unchanged): adjacent-channel
    # MAD / sqrt(2) over every channel of the crop; line slopes raise it, channel correlation lowers it.
    # "line_free": per-pixel robust std of the voxels without line emission (inputs.jl line_free_noise_map);
    # ROHSA gets its median unless noise-per-pixel = true, then the map. noise-channels ("1-10,95-106") replaces the
    # emission mask with fixed line-free channels. A noise-map file wins over both.
    "noise-estimate" => "channel_difference",
    "noise-per-pixel" => "false",
    "noise-channels" => "",
    # Floor of the per-pixel line-free map at this quantile of its values over the pixels with data (0 = no floor).
    # In mosaics, some pixels hold much quieter data in their line-free channels (e.g. single-dish only at the band
    # edges) than in their line channels; the line-free estimate is then too low there and ROHSA would give those
    # pixels' line channels far too much weight. Used only with noise-estimate = "line_free".
    "noise-floor-quantile" => "0",
    # How missing data reach ROHSA. "zero" (default, unchanged): NaN voxels are written as 0 and fitted like data, so
    # a pixel beyond the survey edge is fitted as a zero spectrum. "skip": each pixel's live channel range (first to
    # last channel that is neither 0 nor NaN, from the raw cube) goes to ROHSA in live_channel_range.dat
    # (`filename_range`); channels outside it carry no data term, and pixels without data get noise 0 (no data term;
    # their parameters follow the smoothing term). Needs the cube in memory (write-region-fits = true).
    "missing-data" => "zero",
    # With missing_data = "skip": FITS of reliable channel ranges (rohsa/cli/channel_ranges.jl, HDUs FIRST and LAST,
    # crop-sized) that narrows each pixel's live range, e.g. to drop end-of-band channels a mosaic attenuated. Empty =
    # the live range (first to last channel neither 0 nor NaN).
    "channel-range-file" => "",
    "auto-config-enabled" => "false",
    "auto-config-profile" => "physical_decomposition",
    "auto-config-max-candidates" => "24",
    "auto-config-small-component-fraction" => "0.01",
    "auto-config-broad-fraction-limit" => "0.60",
    "auto-config-sigma-bound-tolerance" => "0.05",
    "auto-config-duplicate-velocity-kms" => "2.0",
    "auto-config-duplicate-sigma-kms" => "1.0",
    "auto-config-physical-weight" => "1.0",
    "auto-config-stability-weight" => "0.05",
)

const REQUIRED_ROHSA_OPTIONS = setdiff(Set(values(ROHSA_CONFIG_KEYS)), Set(keys(OPTIONAL_ROHSA_OPTIONS)))

function config_value_string(value)
    if value isa AbstractVector
        return join(string.(value), ",")
    end
    return string(value)
end

function parse_bool(raw)
    value = lowercase(string(raw))
    value in ("true", "yes", "1", ".true.", "t") && return true
    value in ("false", "no", "0", ".false.", "f") && return false
    throw(ArgumentError("Expected boolean string, got $(repr(raw))"))
end

fortran_bool(flag::Bool) = flag ? ".true." : ".false."

list_items(raw::AbstractString) = [strip(part) for part in split(raw, ",") if !isempty(strip(part))]

# Evaluation counts for checkpoint_at: at most 64 (mod_array.f90), positive, written in ascending order.
function parse_checkpoint_at(raw::AbstractString)
    vals = parse.(Int, list_items(raw))
    length(vals) <= 64 || throw(ArgumentError("checkpoint_at takes at most 64 values"))
    all(>(0), vals) || throw(ArgumentError("checkpoint_at values must be positive"))
    return sort(unique(vals))
end

function parse_int_list(raw::AbstractString)
    vals = parse.(Int, list_items(raw))
    isempty(vals) && throw(ArgumentError("Expected at least one integer in list, got $(repr(raw))"))
    return vals
end

# Settings of the optional `[auto_config]` candidate search. Each axis lists the
# values to try; an axis the config leaves out holds the single `[rohsa]` value.
struct AutoConfig
    enabled::Bool
    profile::String
    max_candidates::Int
    n_gauss::Vector{Int}
    lambda_amp::Vector{Float64}
    lambda_mu::Vector{Float64}
    lambda_sig::Vector{Float64}
    lambda_var_amp::Vector{Float64}
    lambda_var_mu::Vector{Float64}
    lambda_var_sig::Vector{Float64}
    small_component_fraction::Float64
    broad_fraction_limit::Float64
    sigma_bound_tolerance::Float64
    duplicate_velocity_kms::Float64
    duplicate_sigma_kms::Float64
    physical_weight::Float64
    stability_weight::Float64
end

# A ROHSA config, parsed once from the TOML values by `RohsaConfig(opts)`. Code
# reads these fields and never parses option strings again.
Base.@kwdef struct RohsaConfig
    input::String
    output_dir::String
    crop_size::Int
    crop_x0::Int
    crop_y0::Int
    spatial_stride::Int
    reuse::Bool
    write_visuals::Bool
    derive_products::Bool
    write_region_fits::Bool
    write_cubes::Bool
    keep_rohsa_input::Bool
    preflight_only::Bool
    force_full_cube::Bool
    memory_gb::Union{Nothing, Float64}  # nothing: the machine's RAM
    safe_memory_fraction::Float64
    fits_extension::String
    spectral_kernel::String
    clip_negative::Bool
    noise_map::String
    holdout_mask::String
    noise_estimate::String
    noise_per_pixel::Bool
    noise_channels::String
    noise_floor_quantile::Float64
    missing_data::String
    channel_range_file::String
    n_gauss::Vector{Int}
    lambda_amp::Float64
    lambda_mu::Float64
    lambda_sig::Float64
    lambda_var_amp::Float64
    lambda_var_mu::Float64
    lambda_var_sig::Float64
    amp_fact_init::Float64
    sig_init::Float64
    lb_sig_init::Float64
    ub_sig_init::Float64
    lb_sig::Float64
    ub_sig::Float64
    regul_grad_adjoint::Bool
    init_option::String
    maxiter_init::Int
    maxiter::Int
    maxiter_last::Int  # 0: the full-resolution level keeps maxiter
    checkpoint_at::Vector{Int}
    m::Int
    noise::Bool
    regul::Bool
    descent::Bool
    lstd::Int
    ustd::Int
    iprint::Int
    iprint_init::Int
    save_grid::Bool
    auto::AutoConfig
    config_path::String
    # The flattened TOML values as written, recorded verbatim in manifest.toml
    # (a later run reads its geometry keys back from there). No code reads
    # values from here.
    manifest_options::Dict{String, String}
end

function RohsaConfig(o::Dict{String, String})
    getint(key) = parse(Int, o[key])
    getfloat(key) = parse(Float64, o[key])
    getflag(key) = parse_bool(o[key])
    axis(key, base) = haskey(o, "auto-config-$key") ? o["auto-config-$key"] : base
    float_axis(key) = parse.(Float64, list_items(axis(key, o[key])))
    n_gauss = parse_int_list(o["n-gauss"])
    auto = AutoConfig(
        getflag("auto-config-enabled"),
        o["auto-config-profile"],
        getint("auto-config-max-candidates"),
        parse.(Int, list_items(axis("n-gauss", o["n-gauss"]))),
        float_axis("lambda-amp"),
        float_axis("lambda-mu"),
        float_axis("lambda-sig"),
        float_axis("lambda-var-amp"),
        float_axis("lambda-var-mu"),
        float_axis("lambda-var-sig"),
        getfloat("auto-config-small-component-fraction"),
        getfloat("auto-config-broad-fraction-limit"),
        getfloat("auto-config-sigma-bound-tolerance"),
        getfloat("auto-config-duplicate-velocity-kms"),
        getfloat("auto-config-duplicate-sigma-kms"),
        getfloat("auto-config-physical-weight"),
        getfloat("auto-config-stability-weight"),
    )
    return RohsaConfig(;
        input=o["input"],
        output_dir=o["output-dir"],
        crop_size=getint("crop-size"),
        crop_x0=getint("crop-x0"),
        crop_y0=getint("crop-y0"),
        spatial_stride=getint("spatial-stride"),
        reuse=getflag("reuse"),
        write_visuals=getflag("write-visuals"),
        derive_products=getflag("derive-products"),
        write_region_fits=getflag("write-region-fits"),
        write_cubes=getflag("write-cubes"),
        keep_rohsa_input=getflag("keep-rohsa-input"),
        preflight_only=getflag("preflight-only"),
        force_full_cube=getflag("force-full-cube"),
        memory_gb=lowercase(o["memory-gb"]) == "auto" ? nothing : getfloat("memory-gb"),
        safe_memory_fraction=getfloat("safe-memory-fraction"),
        fits_extension=o["fits-extension"],
        spectral_kernel=o["spectral-kernel"],
        clip_negative=getflag("clip-negative"),
        noise_map=o["noise-map"],
        holdout_mask=o["holdout-mask"],
        noise_estimate=o["noise-estimate"],
        noise_per_pixel=getflag("noise-per-pixel"),
        noise_channels=o["noise-channels"],
        noise_floor_quantile=getfloat("noise-floor-quantile"),
        missing_data=o["missing-data"],
        channel_range_file=o["channel-range-file"],
        n_gauss,
        lambda_amp=getfloat("lambda-amp"),
        lambda_mu=getfloat("lambda-mu"),
        lambda_sig=getfloat("lambda-sig"),
        lambda_var_amp=getfloat("lambda-var-amp"),
        lambda_var_mu=getfloat("lambda-var-mu"),
        lambda_var_sig=getfloat("lambda-var-sig"),
        amp_fact_init=getfloat("amp-fact-init"),
        sig_init=getfloat("sig-init"),
        lb_sig_init=getfloat("lb-sig-init"),
        ub_sig_init=getfloat("ub-sig-init"),
        lb_sig=getfloat("lb-sig"),
        ub_sig=getfloat("ub-sig"),
        regul_grad_adjoint=getflag("regul-grad-adjoint"),
        init_option=o["init-option"],
        maxiter_init=getint("maxiter-init"),
        maxiter=getint("maxiter"),
        maxiter_last=getint("maxiter-last"),
        checkpoint_at=parse_checkpoint_at(o["checkpoint-at"]),
        m=getint("m"),
        noise=getflag("rohsa-noise"),
        regul=getflag("regul"),
        descent=getflag("descent"),
        lstd=getint("lstd"),
        ustd=getint("ustd"),
        iprint=getint("iprint"),
        iprint_init=getint("iprint-init"),
        save_grid=getflag("save-grid"),
        auto,
        config_path=o["config-path"],
        manifest_options=o,
    )
end

function load_rohsa_config(path::AbstractString)
    isfile(path) || throw(ArgumentError("ROHSA config file not found: $path"))
    raw = TOML.parsefile(path)
    opts = Dict{String, String}()
    # Flatten TOML tables into the option names used by the CLI and by the
    # generated Fortran namelist, while rejecting silent misspellings.
    for (table_name, table_value) in raw
        table_value isa AbstractDict ||
            throw(ArgumentError("Top-level ROHSA config entry `$table_name` must be a TOML table"))
        for (key, value) in table_value
            if table_name == "auto_config"
                haskey(ROHSA_AUTO_CONFIG_KEYS, key) ||
                    throw(ArgumentError("Unknown ROHSA auto-config key `$table_name.$key` in $path"))
                opts[ROHSA_AUTO_CONFIG_KEYS[key]] = config_value_string(value)
            else
                haskey(ROHSA_CONFIG_KEYS, key) ||
                    throw(ArgumentError("Unknown ROHSA config key `$table_name.$key` in $path"))
                opts[ROHSA_CONFIG_KEYS[key]] = config_value_string(value)
            end
        end
    end
    missing = sort([key for key in REQUIRED_ROHSA_OPTIONS if !haskey(opts, key)])
    isempty(missing) ||
        throw(ArgumentError("ROHSA config $path is missing required options: $(join(missing, ", "))"))
    for (key, value) in OPTIONAL_ROHSA_OPTIONS
        get!(opts, key, value)
    end
    opts["fits-extension"] = VGT.normalize_fits_extension(opts["fits-extension"])
    opts["config-path"] = abspath(path)
    return RohsaConfig(opts)
end

# The science parameters that vary between the runs of one config: the Gaussian
# count and the six regularisation weights (the auto-config axes).
struct RohsaCandidate
    id::String  # empty outside auto-config
    n_gauss::Int
    lambda_amp::Float64
    lambda_mu::Float64
    lambda_sig::Float64
    lambda_var_amp::Float64
    lambda_var_mu::Float64
    lambda_var_sig::Float64
end

function write_rohsa_parameters(path, input_path, rms_path, out_path, timeout_path, cfg::RohsaConfig, cand::RohsaCandidate;
                                range_path::AbstractString="")
    open(path, "w") do io
        # Upstream ROHSA consumes a Fortran namelist. Numeric real values are
        # emitted with d0 exponents to match the original examples.
        println(io, "&user_parameters")
        println(io, "    filename = '$(abspath(input_path))'")
        println(io, "    ,fileout = '$(abspath(out_path))'")
        println(io, "    ,timeout = '$(abspath(timeout_path))'")
        println(io, "    ,filename_noise = '$(abspath(rms_path))'")
        println(io, "    ,n_gauss = $(cand.n_gauss)")
        println(io, "    ,lambda_amp = $(cand.lambda_amp)d0")
        println(io, "    ,lambda_mu = $(cand.lambda_mu)d0")
        println(io, "    ,lambda_sig = $(cand.lambda_sig)d0")
        println(io, "    ,lambda_var_amp = $(cand.lambda_var_amp)d0")
        println(io, "    ,lambda_var_mu = $(cand.lambda_var_mu)d0")
        println(io, "    ,lambda_var_sig = $(cand.lambda_var_sig)d0")
        println(io, "    ,amp_fact_init = $(cfg.amp_fact_init)d0")
        println(io, "    ,sig_init = $(cfg.sig_init)d0")
        println(io, "    ,lb_sig_init = $(cfg.lb_sig_init)d0")
        println(io, "    ,ub_sig_init = $(cfg.ub_sig_init)d0")
        println(io, "    ,lb_sig = $(cfg.lb_sig)d0")
        println(io, "    ,ub_sig = $(cfg.ub_sig)d0")
        println(io, "    ,init_option = '$(cfg.init_option)'")
        println(io, "    ,maxiter_init = $(cfg.maxiter_init)")
        println(io, "    ,maxiter = $(cfg.maxiter)")
        println(io, "    ,m = $(cfg.m)")
        println(io, "    ,noise = $(fortran_bool(cfg.noise))")
        println(io, "    ,regul = $(fortran_bool(cfg.regul))")
        println(io, "    ,descent = $(fortran_bool(cfg.descent))")
        println(io, "    ,lstd = $(cfg.lstd)")
        println(io, "    ,ustd = $(cfg.ustd)")
        println(io, "    ,iprint = $(cfg.iprint)")
        println(io, "    ,iprint_init = $(cfg.iprint_init)")
        println(io, "    ,save_grid = $(fortran_bool(cfg.save_grid))")
        if cfg.regul_grad_adjoint
            println(io, "    ,regul_grad_adjoint = .true.")
        end
        # Live channel range file (missing_data = "skip"); emitted only when set, so older binaries still read the namelist.
        isempty(range_path) || println(io, "    ,filename_range = '$(abspath(range_path))'")
        cfg.maxiter_last > 0 && println(io, "    ,maxiter_last = $(cfg.maxiter_last)")
        if !isempty(cfg.checkpoint_at)
            ckpt_dir = joinpath(dirname(abspath(out_path)), "checkpoints")
            mkpath(ckpt_dir)
            println(io, "    ,checkpoint_at = $(join(cfg.checkpoint_at, ","))")
            println(io, "    ,checkpoint_prefix = '$(joinpath(ckpt_dir, "$(rohsa_run_label(cand))"))'")
        end
        println(io, "    /")
    end
    return path
end

function rohsa_trace_config(cfg::RohsaConfig, cand::RohsaCandidate)
    return Dict(
        "config_path" => cfg.config_path,
        "default_source" => "rohsa/upstream/ROHSA/ROHSApy.ipynb",
        "n_gauss" => cand.n_gauss,
        "lambda_amp" => cand.lambda_amp,
        "lambda_mu" => cand.lambda_mu,
        "lambda_sig" => cand.lambda_sig,
        "lambda_var_amp" => cand.lambda_var_amp,
        "lambda_var_mu" => cand.lambda_var_mu,
        "lambda_var_sig" => cand.lambda_var_sig,
        "amp_fact_init" => cfg.amp_fact_init,
        "sig_init" => cfg.sig_init,
        "lb_sig_init" => cfg.lb_sig_init,
        "ub_sig_init" => cfg.ub_sig_init,
        "lb_sig" => cfg.lb_sig,
        "ub_sig" => cfg.ub_sig,
        "init_option" => cfg.init_option,
        "maxiter_init" => cfg.maxiter_init,
        "maxiter" => cfg.maxiter,
        "m" => cfg.m,
        "noise" => fortran_bool(cfg.noise),
        "regul" => fortran_bool(cfg.regul),
        "descent" => fortran_bool(cfg.descent),
        "lstd" => cfg.lstd,
        "ustd" => cfg.ustd,
        "save_grid" => fortran_bool(cfg.save_grid),
    )
end
