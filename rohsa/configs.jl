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
    # Exact gradient of ROHSA's smoothing term at map edges (mod_array.f90 regul_grad_adjoint). Off by
    # default so existing runs stay reproducible; the original code is wrong in the two outer rows and
    # columns of every level. Emitted only when true, so older ROHSA binaries still read the namelist.
    "regul-grad-adjoint" => "false",
    # Per-pixel noise for ROHSA (2D FITS, crop-sized, Julia x-y order). Empty = one level from the cube
    # (robust_channel_difference_noise). Pixels with NaN or <= 0 are written as 0, which ROHSA skips in the data
    # term at the full-resolution level (their parameters follow the smoothing term): used for held-out blocks.
    "noise-map" => "",
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
    opts["fits-extension"] = Astronomy.normalize_fits_extension(opts["fits-extension"])
    opts["config-path"] = abspath(path)
    return opts
end

function fortran_bool(raw)
    value = lowercase(string(raw))
    value in ("true", "yes", "1", ".true.", "t") && return ".true."
    value in ("false", "no", "0", ".false.", "f") && return ".false."
    throw(ArgumentError("Expected boolean-like value, got $(repr(raw))"))
end

function parse_bool(raw)
    value = lowercase(string(raw))
    value in ("true", "yes", "1") && return true
    value in ("false", "no", "0") && return false
    throw(ArgumentError("Expected boolean string, got $(repr(raw))"))
end

function write_rohsa_parameters(path, input_path, rms_path, out_path, timeout_path, n_gauss::Int, opts)
    open(path, "w") do io
        # Upstream ROHSA consumes a Fortran namelist. Numeric real values are
        # emitted with d0 exponents to match the original examples.
        println(io, "&user_parameters")
        println(io, "    filename = '$(abspath(input_path))'")
        println(io, "    ,fileout = '$(abspath(out_path))'")
        println(io, "    ,timeout = '$(abspath(timeout_path))'")
        println(io, "    ,filename_noise = '$(abspath(rms_path))'")
        println(io, "    ,n_gauss = $n_gauss")
        println(io, "    ,lambda_amp = $(parse(Float64, opts["lambda-amp"]))d0")
        println(io, "    ,lambda_mu = $(parse(Float64, opts["lambda-mu"]))d0")
        println(io, "    ,lambda_sig = $(parse(Float64, opts["lambda-sig"]))d0")
        println(io, "    ,lambda_var_amp = $(parse(Float64, opts["lambda-var-amp"]))d0")
        println(io, "    ,lambda_var_mu = $(parse(Float64, opts["lambda-var-mu"]))d0")
        println(io, "    ,lambda_var_sig = $(parse(Float64, opts["lambda-var-sig"]))d0")
        println(io, "    ,amp_fact_init = $(parse(Float64, opts["amp-fact-init"]))d0")
        println(io, "    ,sig_init = $(parse(Float64, opts["sig-init"]))d0")
        println(io, "    ,lb_sig_init = $(parse(Float64, opts["lb-sig-init"]))d0")
        println(io, "    ,ub_sig_init = $(parse(Float64, opts["ub-sig-init"]))d0")
        println(io, "    ,lb_sig = $(parse(Float64, opts["lb-sig"]))d0")
        println(io, "    ,ub_sig = $(parse(Float64, opts["ub-sig"]))d0")
        println(io, "    ,init_option = '$(opts["init-option"])'")
        println(io, "    ,maxiter_init = $(parse(Int, opts["maxiter-init"]))")
        println(io, "    ,maxiter = $(parse(Int, opts["maxiter"]))")
        println(io, "    ,m = $(parse(Int, opts["m"]))")
        println(io, "    ,noise = $(fortran_bool(opts["rohsa-noise"]))")
        println(io, "    ,regul = $(fortran_bool(opts["regul"]))")
        println(io, "    ,descent = $(fortran_bool(opts["descent"]))")
        println(io, "    ,lstd = $(parse(Int, opts["lstd"]))")
        println(io, "    ,ustd = $(parse(Int, opts["ustd"]))")
        println(io, "    ,iprint = $(parse(Int, opts["iprint"]))")
        println(io, "    ,iprint_init = $(parse(Int, opts["iprint-init"]))")
        println(io, "    ,save_grid = $(fortran_bool(opts["save-grid"]))")
        if parse_bool(get(opts, "regul-grad-adjoint", "false"))
            println(io, "    ,regul_grad_adjoint = .true.")
        end
        println(io, "    /")
    end
    return path
end

function rohsa_trace_config(opts, n_gauss::Int)
    return Dict(
        "config_path" => get(opts, "config-path", ""),
        "default_source" => "rohsa/upstream/ROHSA/ROHSApy.ipynb",
        "n_gauss" => n_gauss,
        "lambda_amp" => parse(Float64, opts["lambda-amp"]),
        "lambda_mu" => parse(Float64, opts["lambda-mu"]),
        "lambda_sig" => parse(Float64, opts["lambda-sig"]),
        "lambda_var_amp" => parse(Float64, opts["lambda-var-amp"]),
        "lambda_var_mu" => parse(Float64, opts["lambda-var-mu"]),
        "lambda_var_sig" => parse(Float64, opts["lambda-var-sig"]),
        "amp_fact_init" => parse(Float64, opts["amp-fact-init"]),
        "sig_init" => parse(Float64, opts["sig-init"]),
        "lb_sig_init" => parse(Float64, opts["lb-sig-init"]),
        "ub_sig_init" => parse(Float64, opts["ub-sig-init"]),
        "lb_sig" => parse(Float64, opts["lb-sig"]),
        "ub_sig" => parse(Float64, opts["ub-sig"]),
        "init_option" => opts["init-option"],
        "maxiter_init" => parse(Int, opts["maxiter-init"]),
        "maxiter" => parse(Int, opts["maxiter"]),
        "m" => parse(Int, opts["m"]),
        "noise" => fortran_bool(opts["rohsa-noise"]),
        "regul" => fortran_bool(opts["regul"]),
        "descent" => fortran_bool(opts["descent"]),
        "lstd" => parse(Int, opts["lstd"]),
        "ustd" => parse(Int, opts["ustd"]),
        "save_grid" => fortran_bool(opts["save-grid"]),
    )
end
