using DrWatson
Base.exit_on_sigint(false)
quickactivate(@__DIR__, "Astronomy")

using CairoMakie
using Dates
using FITSIO
using Printf
using TOML
using Statistics
using Astronomy

const DEFAULT_CONFIG_PATH = normpath(joinpath(@__DIR__, "..", "config.toml"))

Base.include(@__MODULE__, joinpath(@__DIR__, "..", "configs.jl"))
Base.include(@__MODULE__, joinpath(@__DIR__, "..", "memory.jl"))
Base.include(@__MODULE__, joinpath(@__DIR__, "..", "inputs.jl"))
Base.include(@__MODULE__, joinpath(@__DIR__, "..", "runner.jl"))
Base.include(@__MODULE__, joinpath(@__DIR__, "..", "selection.jl"))
Base.include(@__MODULE__, joinpath(@__DIR__, "..", "products_stage.jl"))
Base.include(@__MODULE__, joinpath(@__DIR__, "..", "visuals.jl"))
Base.include(@__MODULE__, joinpath(@__DIR__, "..", "report.jl"))

function parse_args(args)
    config_path = DEFAULT_CONFIG_PATH
    if isempty(args)
        return load_rohsa_config(config_path)
    end
    length(args) == 2 || throw(ArgumentError("Usage: julia --project=. rohsa/cli/preprocess.jl [--config path/to/config.toml]"))
    args[1] == "--config" || throw(ArgumentError("Unknown option $(args[1]). Only --config is accepted; edit rohsa/config.toml for run values."))
    config_path = args[2]
    return load_rohsa_config(config_path)
end

function parse_int_list(raw::AbstractString)
    vals = [parse(Int, strip(part)) for part in split(raw, ",") if !isempty(strip(part))]
    isempty(vals) && throw(ArgumentError("Expected at least one integer in list, got $(repr(raw))"))
    return vals
end

function print_rohsa_preflight(regions, nv::Int, n_gauss_values, opts)
    budget_gib = memory_budget_gib(opts)
    m = parse(Int, opts["m"])
    max_n = maximum(n_gauss_values)
    println("ROHSA memory budget: ", @sprintf("%.1f GiB", budget_gib),
            " (memory_gb=", opts["memory-gb"], ", safe_memory_fraction=", opts["safe-memory-fraction"], ")")
    unsafe = false
    for region in regions
        estimate = rohsa_memory_estimate(region.nx, region.ny, nv, max_n, m)
        @printf("ROHSA preflight [%s]: output=%dx%dx%d padded=%dx%d N=%d m=%d estimated peak %.1f GiB\n",
                region.label, region.nx, region.ny, nv, estimate.padded_side, estimate.padded_side,
                max_n, m, estimate.peak_gib)
        if estimate.peak_gib > budget_gib
            unsafe = true
            suggested = suggested_spatial_stride(region.source_nx, region.source_ny, nv, max_n, m, budget_gib)
            suggested_estimate = rohsa_memory_estimate(cld(region.source_nx, suggested),
                                                       cld(region.source_ny, suggested), nv, max_n, m)
            @printf("Suggested laptop config: spatial_stride = %d, derive_products = false, write_visuals = false (estimated %.1f GiB)\n",
                    suggested, suggested_estimate.peak_gib)
        end
    end
    return unsafe
end

function run_rohsa_preprocessing(opts)
    output_dir = opts["output-dir"]
    mkpath(output_dir)
    input = opts["input"]
    isfile(input) || throw(ArgumentError("Input FITS cube not found: $input"))
    crop_size = parse(Int, opts["crop-size"])
    n_gauss_values = parse_int_list(opts["n-gauss"])
    reuse = parse_bool(opts["reuse"])
    derive_products = parse_bool(opts["derive-products"])
    write_region_fits = parse_bool(opts["write-region-fits"])
    write_visuals = parse_bool(opts["write-visuals"])

    println("ROHSA preprocessing output: ", output_dir)
    println("Input cube: ", input)
    region_label = crop_size <= 0 ? "full cube" : "crop x=$(opts["crop-x0"]) y=$(opts["crop-y0"]) size=$crop_size"
    parse(Int, opts["spatial-stride"]) == 1 || (region_label *= " stride=$(opts["spatial-stride"])")
    println("Region: ", region_label, "  N: ", join(n_gauss_values, ","))
    flush(stdout)

    f = FITS(input)
    header, nx, ny, nv = try
        header = read_header(f[1])
        dims = size(f[1])
        (header, dims...)
    finally
        close(f)
    end
    regions = select_input_regions(input, nx, ny, opts; reuse=reuse)
    println("Selected region: ", only(regions).label)
    flush(stdout)
    unsafe = print_rohsa_preflight(regions, nv, n_gauss_values, opts)
    if parse_bool(opts["preflight-only"])
        println("Preflight-only mode enabled; not writing inputs or running ROHSA.")
        return ""
    end
    if unsafe && !parse_bool(opts["force-full-cube"])
        throw(ArgumentError("ROHSA preflight exceeds the configured memory budget. Use spatial_stride for a reduced-resolution full-footprint run, or force_full_cube = true only on a machine with enough RAM."))
    end
    rohsa_bin = ensure_rohsa_binary()
    println("ROHSA binary: ", rohsa_bin)

    region_rows = Dict{String, Any}[]
    component_rows = Dict{String, Any}[]
    rohsa_metrics = Dict{String, Any}[]
    rohsa_product_rows = Dict{String, Any}[]
    selected_product_rows = Dict{String, Any}[]
    visual_paths = String[]

    for region in regions
        region_dir = joinpath(output_dir, region.label)
        println("[", region.label, "] writing ROHSA input products")
        flush(stdout)
        retain_cube = derive_products || write_visuals
        region_products = write_region_products(input, header, region, region_dir;
                                                retain_cube=retain_cube,
                                                write_region_fits=write_region_fits)
        push!(region_rows, Dict{String, Any}(
            "region" => region.label,
            "x0" => region.x0,
            "y0" => region.y0,
            "nx" => region.nx,
            "ny" => region.ny,
            "spatial_stride" => region.x_stride,
            "noise_estimate" => region_products.sigma,
            "fits_path" => region_products.fits_path,
        ))

        region_rohsa_metrics = Dict{String, Any}[]
        region_product_rows = Dict{String, Any}[]
        for n_gauss in n_gauss_values
            println("[", region.label, "] running ROHSA N=", n_gauss)
            flush(stdout)
            run_stage = run_rohsa_decomposition_stage!(
                region_dir,
                region,
                region_products,
                n_gauss,
                opts,
                rohsa_bin,
                reuse,
            )
            if !derive_products
                push!(region_rohsa_metrics, Dict{String, Any}(
                    "region" => region.label,
                    "n_gauss" => n_gauss,
                    "status" => run_stage.status,
                    "selected_for_products" => false,
                    "run_trace" => run_stage.trace_path,
                    "products_trace" => "",
                    "rohsa_output" => run_stage.rohsa_out,
                    "product_root" => "",
                ))
                continue
            end
            written = derive_rohsa_products_stage(region_dir, region.label, region_products.header,
                                                  region_products.cube, run_stage.rohsa_out, n_gauss)
            for row in written.component_rows
                out = copy(row)
                out["fits_path"] = written.component_paths[Int(row["component"])]
                push!(component_rows, out)
            end
            metric_row = copy(written.fit_row)
            metric_row["selection_score"] = rohsa_selection_score(metric_row)
            metric_row["selected_for_products"] = false
            metric_row["product_root"] = written.root
            metric_row["run_trace"] = run_stage.trace_path
            metric_row["products_trace"] = written.trace_path
            push!(region_rohsa_metrics, metric_row)
            for (label, path) in written.map_paths
                push!(region_product_rows, Dict{String, Any}(
                    "region" => region.label,
                    "n_gauss" => n_gauss,
                    "kind" => "map",
                    "label" => label,
                    "path" => path,
                ))
            end
            for (label, path) in written.cube_paths
                push!(region_product_rows, Dict{String, Any}(
                    "region" => region.label,
                    "n_gauss" => n_gauss,
                    "kind" => "cube",
                    "label" => label,
                    "path" => path,
                ))
            end
            if write_visuals
                println("[", region.label, "] writing ROHSA N=", n_gauss, " visuals")
                flush(stdout)
                append!(visual_paths, write_rohsa_visuals(region_dir, n_gauss,
                                                          region_products.cube, written))
            end
        end
        if derive_products
            best_metric = select_best_rohsa_n(region_rohsa_metrics)
            selected_n = Int(best_metric["n_gauss"])
            for row in region_rohsa_metrics
                row["selected_for_products"] = Int(row["n_gauss"]) == selected_n
            end
            append!(selected_product_rows, [row for row in region_product_rows if Int(row["n_gauss"]) == selected_n])
            @printf("[%s] selected ROHSA N=%d for derived products (score %.4f, product root %s)\n",
                    region.label, selected_n, best_metric["selection_score"], best_metric["product_root"])
        else
            println("[", region.label, "] product derivation disabled; no ROHSA configuration selected")
        end
        append!(rohsa_metrics, region_rohsa_metrics)
        append!(rohsa_product_rows, region_product_rows)
        flush(stdout)
    end

    println("Writing summary tables and report")
    paths = String[]
    selected_rows = Dict{String, Any}[]
    selection_path = ""
    if derive_products
        selected_by_region = select_best_rohsa_by_group(rohsa_metrics; group_key="region")
        selected_rows = selected_rohsa_summary_rows(selected_by_region; group_label="region")
        selection_path = write_rohsa_selection(joinpath(output_dir, "selected_rohsa_config.toml"),
                                               selected_by_region)
    end
    push!(paths, write_csv(joinpath(output_dir, "region_selection.csv"), region_rows,
                           ["region", "x0", "y0", "nx", "ny", "spatial_stride", "noise_estimate", "fits_path"]))
    push!(paths, write_csv(joinpath(output_dir, "rohsa_fit_metrics.csv"), rohsa_metrics,
                           ["region", "n_gauss", "status", "selected_for_products", "selection_score", "total_model_over_data", "voxel_rmse", "voxel_mae", "voxel_residual_skew", "moment_relative_rmse", "moment_mean_relative_bias", "moment_residual_skew", "rohsa_output", "product_root", "run_trace", "products_trace"]))
    push!(paths, write_csv(joinpath(output_dir, "rohsa_component_summary.csv"), component_rows,
                           ["region", "n_gauss", "component", "width_group", "td_group", "inband_emission_fraction", "weighted_sigma_kms", "weighted_mu_kms", "weighted_t_doppler_k", "amp_fill_fraction_gt_5pct_peak", "fits_path"]))
    push!(paths, write_csv(joinpath(output_dir, "rohsa_product_paths.csv"), rohsa_product_rows,
                           ["region", "n_gauss", "kind", "label", "path"]))
    push!(paths, write_csv(joinpath(output_dir, "selected_rohsa_products.csv"), selected_product_rows,
                           ["region", "n_gauss", "kind", "label", "path"]))
    push!(paths, write_csv(joinpath(output_dir, "rohsa_selection_summary.csv"), selected_rows,
                           ["region", "selected_n_gauss", "selection_score", "total_model_over_data", "moment_relative_rmse", "moment_mean_relative_bias", "moment_residual_skew", "product_root", "run_trace", "products_trace"]))
    !isempty(selection_path) && push!(paths, selection_path)
    push!(paths, write_manifest(joinpath(output_dir, "manifest.toml"), opts, rohsa_bin, visual_paths))
    report_target = joinpath(output_dir, "report.md")
    visual_index_path = write_visuals ? write_visual_index(joinpath(output_dir, "visual_index.html"), visual_paths, vcat(paths, [report_target])) : ""
    write_visuals && push!(paths, visual_index_path)
    report_path = write_report(report_target, opts, region_rows, rohsa_metrics, selected_rows, paths, visual_paths)
    println("Wrote primary report: ", report_path)
    println("Selected products manifest: ", joinpath(output_dir, "selected_rohsa_products.csv"))
    println("All ROHSA products saved under: ", output_dir)
    return report_path
end

if abspath(PROGRAM_FILE) == @__FILE__
    try
        run_rohsa_preprocessing(parse_args(ARGS))
    catch err
        showerror(stderr, err)
        println(stderr)
        exit(1)
    end
end
