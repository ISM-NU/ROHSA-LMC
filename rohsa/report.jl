function csv_escape(value)
    text = string(value)
    if occursin(",", text) || occursin("\"", text) || occursin("\n", text)
        return string("\"", replace(text, "\"" => "\"\""), "\"")
    end
    return text
end

function write_csv(path::AbstractString, rows::Vector{<:AbstractDict}, columns)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, join(columns, ","))
        for row in rows
            println(io, join([csv_escape(get(row, col, "")) for col in columns], ","))
        end
    end
    return path
end

function sorted_rows(rows, columns)
    return sort(collect(rows); by=row -> tuple([string(get(row, col, "")) for col in columns]...))
end

function write_report(path, opts, region_rows, rohsa_metrics, selected_rows, output_paths, visual_paths)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "# ROHSA Preprocessing Pipeline")
        println(io)
        println(io, "Generated UTC: `", Dates.format(Dates.now(Dates.UTC), dateformat"yyyy-mm-ddTHH:MM:SS"), "Z`")
        println(io)
        println(io, "This derives ROHSA decomposition products only. It does not run or compare VGT.")
        println(io)
        println(io, "## Inputs")
        println(io)
        println(io, "- H I cube: `", opts["input"], "`")
        println(io, "- ROHSA Gaussian counts: `", opts["n-gauss"], "`")
        println(io)
        println(io, "## Input Region")
        println(io)
        println(io, "| Region | x0 | y0 | nx | ny |")
        println(io, "|---|---:|---:|---:|---:|")
        for row in sorted_rows(region_rows, ["region"])
            @printf(io, "| `%s` | %d | %d | %d | %d |\n",
                    row["region"], row["x0"], row["y0"], row["nx"], row["ny"])
        end
        println(io)
        println(io, "## ROHSA Fit Metrics")
        println(io)
        if parse_bool(get(opts, "derive-products", "true"))
            if parse_bool(get(opts, "auto-config-enabled", "false"))
                println(io, "| Region | Candidate | N | selected | score | recon | physical | stability | model/data | moment RMSE | moment bias | residual skew |")
                println(io, "|---|---|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|")
                for row in sorted_rows(rohsa_metrics, ["region", "candidate_id", "n_gauss"])
                    @printf(io, "| `%s` | `%s` | %d | `%s` | %.4f | %.4f | %.4f | %.4f | %.4f | %.4f | %.4f | %.4f |\n",
                            row["region"], get(row, "candidate_id", ""), row["n_gauss"],
                            row["selected_for_products"], row["selection_score"],
                            get(row, "reconstruction_score", row["selection_score"]),
                            get(row, "physical_penalty", 0.0), get(row, "stability_penalty", 0.0),
                            row["total_model_over_data"], row["moment_relative_rmse"],
                            row["moment_mean_relative_bias"], row["moment_residual_skew"])
                end
            else
                println(io, "| Region | N | selected | score | model/data | moment RMSE | moment bias | residual skew | voxel RMSE |")
                println(io, "|---|---:|---|---:|---:|---:|---:|---:|---:|")
                for row in sorted_rows(rohsa_metrics, ["region", "n_gauss"])
                    @printf(io, "| `%s` | %d | `%s` | %.4f | %.4f | %.4f | %.4f | %.4f | %.3f |\n",
                            row["region"], row["n_gauss"], row["selected_for_products"], row["selection_score"],
                            row["total_model_over_data"], row["moment_relative_rmse"],
                            row["moment_mean_relative_bias"], row["moment_residual_skew"], row["voxel_rmse"])
                end
            end
        else
            println(io, "Product derivation was disabled, so model/residual metrics were not computed.")
            println(io)
            println(io, "| Region | N | ROHSA status | run trace | output |")
            println(io, "|---|---:|---|---|---|")
            for row in sorted_rows(rohsa_metrics, ["region", "n_gauss"])
                @printf(io, "| `%s` | %d | `%s` | `%s` | `%s` |\n",
                        row["region"], row["n_gauss"], row["status"], row["run_trace"], row["rohsa_output"])
            end
        end
        println(io)
        if parse_bool(get(opts, "derive-products", "true"))
            if parse_bool(get(opts, "auto-config-enabled", "false"))
                println(io, "Auto-config selected the lowest physical-decomposition score: reconstruction quality plus penalties for negligible components, sigma-bound pileup, excessive broad-component fraction, duplicate components, and rough parameter maps. This is a preprocessing diagnostic, not downstream VGT validation.")
            else
                println(io, "The selected ROHSA run minimizes a reporting-only score based on flux recovery, moment residuals, moment bias, and residual skewness. Residual skewness near zero is preferred because positive skew suggests missed emission and negative skew suggests overfitting.")
            end
        end
        println(io)
        println(io, "## Selected Products")
        println(io)
        if isempty(selected_rows)
            println(io, "No selected product set was written because product derivation was disabled.")
        else
            if parse_bool(get(opts, "auto-config-enabled", "false"))
                println(io, "| Region | selected candidate | selected N | score | recon | physical | stability | product root |")
                println(io, "|---|---|---:|---:|---:|---:|---:|---|")
            else
                println(io, "| Region | selected N | score | product root |")
                println(io, "|---|---:|---:|---|")
            end
            for row in sorted_rows(selected_rows, ["region"])
                if parse_bool(get(opts, "auto-config-enabled", "false"))
                    @printf(io, "| `%s` | `%s` | %d | %.4f | %.4f | %.4f | %.4f | `%s` |\n",
                            row["region"], get(row, "selected_candidate_id", ""),
                            row["selected_n_gauss"], row["selection_score"],
                            get(row, "reconstruction_score", row["selection_score"]),
                            get(row, "physical_penalty", 0.0),
                            get(row, "stability_penalty", 0.0), row["product_root"])
                else
                    @printf(io, "| `%s` | %d | %.4f | `%s` |\n",
                            row["region"], row["selected_n_gauss"], row["selection_score"], row["product_root"])
                end
            end
        end
        println(io)
        println(io, "## Visual Outputs")
        println(io)
        if isempty(visual_paths)
            println(io, "Visual output generation was disabled.")
        else
            report_dir = dirname(path)
            for visual_path in visual_paths
                rel = relpath(visual_path, report_dir)
                println(io, "![", basename(visual_path), "](", rel, ")")
                println(io)
                println(io, "`", visual_path, "`")
                println(io)
            end
        end
        println(io)
        println(io, "## Output Files")
        println(io)
        for out in output_paths
            println(io, "- `", out, "`")
        end
        println(io)
        println(io, "## Interpretation Rule")
        println(io)
        println(io, "Treat ROHSA outputs as preprocessing products. Any VGT or starlight validation should be run as a separate downstream analysis using the generated product manifest.")
    end
    return path
end

function write_manifest(path::AbstractString, opts, rohsa_bin::AbstractString, visual_paths)
    mkpath(dirname(path))
    open(path, "w") do io
        TOML.print(io, Dict(
            "generated_at_utc" => string(Dates.format(Dates.now(Dates.UTC), dateformat"yyyy-mm-ddTHH:MM:SS"), "Z"),
            "purpose" => "ROHSA preprocessing and product derivation for downstream VGT input cubes; no VGT comparison is run.",
            "options" => opts,
            "rohsa_bin" => rohsa_bin,
            "visual_outputs" => visual_paths,
        ); sorted=true)
    end
    return path
end
