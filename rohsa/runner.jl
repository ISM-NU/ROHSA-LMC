upstream_rohsa_source_dir() = joinpath(@__DIR__, "upstream", "ROHSA")
upstream_rohsa_binary() = joinpath(upstream_rohsa_source_dir(), "src", "ROHSA")

function ensure_rohsa_binary()
    bin = upstream_rohsa_binary()
    isfile(bin) && return abspath(bin)
    source_dir = upstream_rohsa_source_dir()
    isdir(source_dir) ||
        throw(ArgumentError("Bundled upstream ROHSA source not found at $source_dir"))
    src_dir = joinpath(source_dir, "src")
    if Sys.which("gfortran") !== nothing
        run(Cmd(`make clean all`, dir=src_dir))
    else
        throw(ArgumentError("No gfortran found. Install gfortran or build $(upstream_rohsa_binary()) manually."))
    end
    isfile(bin) || throw(ErrorException("ROHSA build completed but binary was not found at $bin"))
    return abspath(bin)
end

function run_rohsa(rohsa_bin::String, params_path::String, out_path::String, reuse::Bool)
    if reuse && isfile(out_path)
        return "reused"
    end
    run(`$rohsa_bin $params_path`)
    isfile(out_path) || throw(ErrorException("ROHSA did not write expected output $out_path"))
    return "ran"
end

function run_rohsa_decomposition_stage!(region_dir, region, region_products, n_gauss::Int, opts, rohsa_bin::String, reuse::Bool)
    params_path = joinpath(region_dir, "parameters_n$(n_gauss).txt")
    rohsa_out = joinpath(region_dir, "rohsa_n$(n_gauss).dat")
    timeout_path = joinpath(region_dir, "timestep_n$(n_gauss).dat")
    write_rohsa_parameters(params_path, region_products.dat_path, region_products.rms_path,
                           rohsa_out, timeout_path, n_gauss, opts)
    status = run_rohsa(rohsa_bin, params_path, rohsa_out, reuse)
    trace_path = joinpath(region_dir, "rohsa_n$(n_gauss)_run_trace.toml")
    open(trace_path, "w") do io
        TOML.print(io, Dict(
            "stage" => "rohsa_decomposition",
            "region" => region.label,
            "status" => status,
            "input_dat" => region_products.dat_path,
            "rms_map" => region_products.rms_path,
            "parameters" => params_path,
            "rohsa_output" => rohsa_out,
            "timeout" => timeout_path,
            "config" => rohsa_trace_config(opts, n_gauss),
        ); sorted=true)
    end
    return (; params_path, rohsa_out, timeout_path, trace_path, status)
end
