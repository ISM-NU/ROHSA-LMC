function metric_float(row, key::AbstractString)
    value = row[key]
    value isa Number && return Float64(value)
    return parse(Float64, string(value))
end

function rohsa_selection_score(row)
    # Reporting-only heuristic for choosing which N gets product manifests. It
    # favors flux recovery and low moment residuals; it is not a VGT validation.
    vals = (
        abs(metric_float(row, "total_model_over_data") - 1.0),
        metric_float(row, "moment_relative_rmse"),
        abs(metric_float(row, "moment_mean_relative_bias")),
        0.10 * abs(metric_float(row, "moment_residual_skew")),
    )
    any(!isfinite, vals) && return Inf
    return sum(vals)
end

function select_best_rohsa_n(rows)
    isempty(rows) && throw(ArgumentError("Cannot select best ROHSA run from empty metrics"))
    scores = [haskey(row, "selection_score") ? metric_float(row, "selection_score") : rohsa_selection_score(row) for row in rows]
    return rows[argmin(scores)]
end

function grouped_rohsa_rows(rows; group_key="region")
    grouped = Dict{String, Vector{Dict{String, Any}}}()
    for row in rows
        key = string(get(row, group_key, ""))
        push!(get!(grouped, key, Dict{String, Any}[]), row)
    end
    return grouped
end

function select_best_rohsa_by_group(rows; group_key="region")
    out = Dict{String, Dict{String, Any}}()
    for (key, group_rows) in grouped_rohsa_rows(rows; group_key=group_key)
        best = copy(select_best_rohsa_n(group_rows))
        best["selection_score"] = get(best, "selection_score", rohsa_selection_score(best))
        out[key] = best
    end
    return out
end

function selected_rohsa_summary_rows(selected; group_label="region")
    rows = Dict{String, Any}[]
    for key in sort(collect(keys(selected)))
        row = selected[key]
        push!(rows, Dict{String, Any}(
            group_label => key,
            "selected_n_gauss" => row["n_gauss"],
            "selected_candidate_id" => get(row, "candidate_id", ""),
            "selection_score" => row["selection_score"],
            "reconstruction_score" => get(row, "reconstruction_score", row["selection_score"]),
            "physical_penalty" => get(row, "physical_penalty", 0.0),
            "stability_penalty" => get(row, "stability_penalty", 0.0),
            "total_model_over_data" => row["total_model_over_data"],
            "moment_relative_rmse" => row["moment_relative_rmse"],
            "moment_mean_relative_bias" => row["moment_mean_relative_bias"],
            "moment_residual_skew" => row["moment_residual_skew"],
            "product_root" => row["product_root"],
            "run_trace" => row["run_trace"],
            "products_trace" => row["products_trace"],
        ))
    end
    return rows
end

function write_rohsa_selection(path::AbstractString, selected)
    mkpath(dirname(path))
    open(path, "w") do io
        TOML.print(io, Dict(
            "stage" => "rohsa_config_selection",
            "selection_rule" => "minimize selection_score; auto-config uses reconstruction_score+physical_penalty+stability_penalty",
            "selected" => selected,
        ); sorted=true)
    end
    return path
end
