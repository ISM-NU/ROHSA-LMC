const AUTO_CONFIG_AXIS_KEYS = (
    "n-gauss",
    "lambda-amp",
    "lambda-mu",
    "lambda-sig",
    "lambda-var-amp",
    "lambda-var-mu",
    "lambda-var-sig",
)

const AUTO_CONFIG_COLUMNS = [
    "region",
    "candidate_id",
    "n_gauss",
    "lambda_amp",
    "lambda_mu",
    "lambda_sig",
    "lambda_var_amp",
    "lambda_var_mu",
    "lambda_var_sig",
    "selected_for_products",
    "selection_score",
    "reconstruction_score",
    "physical_penalty",
    "stability_penalty",
    "small_component_penalty",
    "sigma_bound_penalty",
    "broad_fraction_penalty",
    "duplicate_component_penalty",
    "roughness_penalty",
    "total_model_over_data",
    "moment_relative_rmse",
    "moment_mean_relative_bias",
    "moment_residual_skew",
    "rohsa_output",
    "product_root",
    "run_trace",
    "products_trace",
]

function auto_candidate_id(index::Int, n_gauss::Int)
    return @sprintf("c%03d_n%d", index, n_gauss)
end

function limit_auto_candidates(candidates, max_candidates::Int)
    max_candidates >= 1 || throw(ArgumentError("auto_config.max_candidates must be >= 1"))
    length(candidates) <= max_candidates && return candidates
    total = length(candidates)
    selected = Int[]
    for raw_idx in round.(Int, range(1, total; length=max_candidates))
        raw_idx = clamp(raw_idx, 1, total)
        raw_idx in selected || push!(selected, raw_idx)
    end
    idx = 1
    while length(selected) < max_candidates && idx <= total
        idx in selected || push!(selected, idx)
        idx += 1
    end
    sort!(selected)
    return candidates[selected]
end

function build_auto_config_candidates(cfg::RohsaConfig)
    a = cfg.auto
    axes = (a.n_gauss, a.lambda_amp, a.lambda_mu, a.lambda_sig, a.lambda_var_amp, a.lambda_var_mu, a.lambda_var_sig)
    for (key, axis) in zip(AUTO_CONFIG_AXIS_KEYS, axes)
        isempty(axis) && throw(ArgumentError("Auto-config axis `$key` resolved to an empty value"))
    end
    # The last axis varies fastest.
    grid = [reverse(combo) for combo in vec(collect(Iterators.product(reverse(axes)...)))]
    limited = limit_auto_candidates(grid, a.max_candidates)
    return [RohsaCandidate(auto_candidate_id(i, combo[1]), combo...) for (i, combo) in enumerate(limited)]
end

function resolved_rohsa_candidates(cfg::RohsaConfig)
    a = cfg.auto
    if !a.enabled
        length(cfg.n_gauss) == 1 ||
            throw(ArgumentError("n_gauss lists $(length(cfg.n_gauss)) values; use [auto_config] to run several, or give one"))
        return [RohsaCandidate("", only(cfg.n_gauss), cfg.lambda_amp, cfg.lambda_mu, cfg.lambda_sig,
                               cfg.lambda_var_amp, cfg.lambda_var_mu, cfg.lambda_var_sig)]
    end
    a.profile == "physical_decomposition" ||
        throw(ArgumentError("Unsupported auto_config.profile=$(repr(a.profile)); expected \"physical_decomposition\""))
    return build_auto_config_candidates(cfg)
end

function reconstruction_score(row)
    vals = (
        abs(metric_float(row, "total_model_over_data") - 1.0),
        metric_float(row, "moment_relative_rmse"),
        abs(metric_float(row, "moment_mean_relative_bias")),
        0.10 * abs(metric_float(row, "moment_residual_skew")),
    )
    any(!isfinite, vals) && return Inf
    return sum(vals)
end

function finite_float_values(arr)
    return [Float64(v) for v in arr if isfinite(v)]
end

function normalized_roughness(map)
    vals = finite_float_values(map)
    length(vals) < 4 && return 0.0
    scale = quantile(vals, 0.90) - quantile(vals, 0.10)
    if !(isfinite(scale) && scale > 0)
        scale = max(abs(mean(vals)), 1.0)
    end
    rough = Float64[]
    rows, cols = size(map)
    if rows < 3 || cols < 3
        return 0.0
    end
    @inbounds for j in 2:(cols - 1), i in 2:(rows - 1)
        center = Float64(map[i, j])
        neigh = (
            Float64(map[i - 1, j]),
            Float64(map[i + 1, j]),
            Float64(map[i, j - 1]),
            Float64(map[i, j + 1]),
        )
        if isfinite(center) && all(isfinite, neigh)
            push!(rough, abs(4center - sum(neigh)) / scale)
        end
    end
    isempty(rough) && return 0.0
    return min(mean(rough), 10.0)
end

function duplicate_component_penalty(component_rows, duplicate_velocity_kms::Float64,
                                     duplicate_sigma_kms::Float64)
    penalty = 0.0
    for i in 1:length(component_rows), j in (i + 1):length(component_rows)
        row_i = component_rows[i]
        row_j = component_rows[j]
        mu_close = abs(metric_float(row_i, "weighted_mu_kms") -
                       metric_float(row_j, "weighted_mu_kms")) <= duplicate_velocity_kms
        sig_close = abs(metric_float(row_i, "weighted_sigma_kms") -
                        metric_float(row_j, "weighted_sigma_kms")) <= duplicate_sigma_kms
        if mu_close && sig_close
            penalty += min(metric_float(row_i, "inband_emission_fraction"),
                           metric_float(row_j, "inband_emission_fraction"))
        end
    end
    return penalty
end

function physical_score_terms(params, products, cfg::RohsaConfig)
    rows = collect(products.component_rows)
    n_gauss = length(rows)
    a = cfg.auto
    min_fraction = a.small_component_fraction
    broad_limit = a.broad_fraction_limit
    bound_tol = a.sigma_bound_tolerance
    duplicate_velocity = a.duplicate_velocity_kms
    duplicate_sigma = a.duplicate_sigma_kms
    physical_weight = a.physical_weight
    stability_weight = a.stability_weight
    lb_sig = cfg.lb_sig
    ub_sig = cfg.ub_sig
    sig_span = max(ub_sig - lb_sig, eps(Float64))
    low_bound = lb_sig + bound_tol * sig_span
    high_bound = ub_sig - bound_tol * sig_span

    fractions = [metric_float(row, "inband_emission_fraction") for row in rows]
    small_component_penalty = isempty(fractions) ? 0.0 :
        sum(max(0.0, min_fraction - fraction) / max(min_fraction, eps(Float64)) for fraction in fractions) / n_gauss

    sigma_bound_hits = 0
    sigma_count = 0
    for g in 1:n_gauss
        sig = @view params[3g, :, :]
        for value in sig
            if isfinite(value)
                sigma_count += 1
                (value <= low_bound || value >= high_bound) && (sigma_bound_hits += 1)
            end
        end
    end
    sigma_bound_penalty = sigma_count == 0 ? 0.0 : sigma_bound_hits / sigma_count

    broad_fraction_vals = finite_float_values(products.group_nhi["broad_fraction"])
    broad_fraction = isempty(broad_fraction_vals) ? 0.0 : mean(broad_fraction_vals)
    broad_fraction_penalty = max(0.0, broad_fraction - broad_limit) / max(1.0 - broad_limit, eps(Float64))

    duplicate_penalty = duplicate_component_penalty(rows, duplicate_velocity, duplicate_sigma)

    weighted_roughness = 0.0
    total_weight = 0.0
    for g in 1:n_gauss
        weight = max(fractions[g], 0.0)
        weight == 0.0 && continue
        amp = @view params[3g - 2, :, :]
        mu = @view params[3g - 1, :, :]
        sig = @view params[3g, :, :]
        component_roughness = mean((
            normalized_roughness(amp),
            normalized_roughness(mu),
            normalized_roughness(sig),
        ))
        weighted_roughness += weight * component_roughness
        total_weight += weight
    end
    roughness_penalty = total_weight == 0.0 ? 0.0 : weighted_roughness / total_weight

    physical_penalty = physical_weight * (
        small_component_penalty +
        sigma_bound_penalty +
        0.5 * broad_fraction_penalty +
        duplicate_penalty
    )
    stability_penalty = stability_weight * roughness_penalty

    return Dict{String, Any}(
        "small_component_penalty" => small_component_penalty,
        "sigma_bound_penalty" => sigma_bound_penalty,
        "broad_fraction_penalty" => broad_fraction_penalty,
        "duplicate_component_penalty" => duplicate_penalty,
        "roughness_penalty" => roughness_penalty,
        "physical_penalty" => physical_penalty,
        "stability_penalty" => stability_penalty,
    )
end

function apply_auto_config_scores!(metric_row, params, products, cfg::RohsaConfig)
    recon = reconstruction_score(metric_row)
    metric_row["reconstruction_score"] = recon
    if cfg.auto.enabled
        physical_terms = physical_score_terms(params, products, cfg)
        for (key, value) in physical_terms
            metric_row[key] = value
        end
        metric_row["selection_score"] = recon + metric_float(metric_row, "physical_penalty") +
                                        metric_float(metric_row, "stability_penalty")
        metric_row["selection_mode"] = "physical_decomposition"
    else
        for key in (
            "small_component_penalty",
            "sigma_bound_penalty",
            "broad_fraction_penalty",
            "duplicate_component_penalty",
            "roughness_penalty",
            "physical_penalty",
            "stability_penalty",
        )
            metric_row[key] = 0.0
        end
        metric_row["selection_score"] = recon
        metric_row["selection_mode"] = "reconstruction"
    end
    return metric_row
end

function candidate_summary_values(cand::RohsaCandidate)
    return Dict(
        "candidate_id" => cand.id,
        "n_gauss" => cand.n_gauss,
        "lambda_amp" => cand.lambda_amp,
        "lambda_mu" => cand.lambda_mu,
        "lambda_sig" => cand.lambda_sig,
        "lambda_var_amp" => cand.lambda_var_amp,
        "lambda_var_mu" => cand.lambda_var_mu,
        "lambda_var_sig" => cand.lambda_var_sig,
    )
end
