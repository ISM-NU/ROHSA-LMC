# ROHSA config loading and parameter-file serialization.
#
# Provenance:
# - `rohsa/config.toml` starts from rohsa/upstream/ROHSA/ROHSApy.ipynb,
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
    "preflight_only" => "preflight-only",
    "force_full_cube" => "force-full-cube",
    "memory_gb" => "memory-gb",
    "safe_memory_fraction" => "safe-memory-fraction",
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

const OPTIONAL_ROHSA_OPTIONS = Dict(
    "spatial-stride" => "1",
    "derive-products" => "true",
    "write-region-fits" => "true",
    "preflight-only" => "false",
    "force-full-cube" => "false",
    "memory-gb" => "auto",
    "safe-memory-fraction" => "0.70",
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
    for (table_name, table_value) in raw
        table_value isa AbstractDict ||
            throw(ArgumentError("Top-level ROHSA config entry `$table_name` must be a TOML table"))
        for (key, value) in table_value
            haskey(ROHSA_CONFIG_KEYS, key) ||
                throw(ArgumentError("Unknown ROHSA config key `$table_name.$key` in $path"))
            opts[ROHSA_CONFIG_KEYS[key]] = config_value_string(value)
        end
    end
    missing = sort([key for key in REQUIRED_ROHSA_OPTIONS if !haskey(opts, key)])
    isempty(missing) ||
        throw(ArgumentError("ROHSA config $path is missing required options: $(join(missing, ", "))"))
    for (key, value) in OPTIONAL_ROHSA_OPTIONS
        get!(opts, key, value)
    end
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

function write_rohsa_parameters(path, dat_path, rms_path, out_path, timeout_path, n_gauss::Int, opts)
    open(path, "w") do io
        println(io, "&user_parameters")
        println(io, "    filename = '$(abspath(dat_path))'")
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
