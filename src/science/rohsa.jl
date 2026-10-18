const HI_COLUMN_DENSITY_FACTOR = 1.82243e18
const DOPPLER_TEMPERATURE_FACTOR = 121.1

function velocity_axis_kms(header, nv::Int)
    crval = Float64(header["CRVAL3"])
    crpix = Float64(header["CRPIX3"])
    cdelt = Float64(header["CDELT3"])
    # FITS convention: channel k (1-based) is at CRVAL3 + (k - CRPIX3) * CDELT3.
    native = [crval + (k - crpix) * cdelt for k in 1:nv]
    unit = try
        lowercase(strip(String(header["CUNIT3"])))
    catch
        ""
    end
    # Some survey cubes omit or vary CUNIT3 spelling; fall back to value scale
    # so native m/s axes are still converted before physical products are made.
    scale = if unit in ("m/s", "ms-1", "m s-1", "m/s ")
        1e-3
    elseif unit in ("km/s", "kms-1", "km s-1")
        1.0
    else
        maximum(abs.(native)) > 1_000 ? 1e-3 : 1.0
    end
    return native .* scale, abs(cdelt * scale)
end

function readdlm_skip_comments(path::AbstractString)
    rows = Vector{Vector{Float64}}()
    for line in eachline(path)
        stripped = strip(line)
        (isempty(stripped) || startswith(stripped, "#")) && continue
        push!(rows, parse.(Float64, split(stripped)))
    end
    isempty(rows) && return zeros(Float64, 0, 0)
    out = zeros(Float64, length(rows), maximum(length, rows))
    for (i, row) in enumerate(rows)
        out[i, 1:length(row)] .= row
    end
    return out
end

function read_rohsa_parameters(path::AbstractString, n_gauss::Int, nx::Int, ny::Int)
    data = readdlm_skip_comments(path)
    expected = nx * ny * n_gauss
    size(data, 1) >= expected ||
        throw(ArgumentError("ROHSA output has $(size(data, 1)) rows; expected at least $expected"))
    params = zeros(Float64, 3n_gauss, nx, ny)
    idx = 1
    # ROHSA writes rows in y, x, Gaussian order with zero-based spectral centers.
    # Store parameters as [amp, mu_channel, sigma_channel] planes per component.
    for y in 1:ny, x in 1:nx, g in 1:n_gauss
        params[3g - 2, x, y] = data[idx, 3]
        params[3g - 1, x, y] = data[idx, 4] - 1.0
        params[3g, x, y] = data[idx, 5]
        idx += 1
    end
    return params
end

function reconstruct_rohsa_components(params, nv::Int)
    n_gauss = size(params, 1) ÷ 3
    nx, ny = size(params, 2), size(params, 3)
    components = [zeros(Float32, nx, ny, nv) for _ in 1:n_gauss]
    for g in 1:n_gauss
        amp = @view params[3g - 2, :, :]
        mu = @view params[3g - 1, :, :]
        sig = @view params[3g, :, :]
        for k in 1:nv
            channel = Float64(k - 1)
            @inbounds for y in 1:ny, x in 1:nx
                # Keep a positive floor to avoid singular fitted widths when
                # reconstructing the Gaussian model cube from parameter maps.
                sigma = max(sig[x, y], 1e-6)
                components[g][x, y, k] = Float32(amp[x, y] * exp(-((channel - mu[x, y])^2) / (2sigma^2)))
            end
        end
    end
    return components
end

moment0(cube) = dropdims(sum(Float64.(cube); dims=3); dims=3)

function weighted_spectral_moment(cube, velocities_kms::AbstractVector)
    nx, ny, nv = size(cube)
    length(velocities_kms) == nv ||
        throw(ArgumentError("Velocity axis length $(length(velocities_kms)) does not match cube channels $nv"))
    out = fill(NaN, nx, ny)
    @inbounds for y in 1:ny, x in 1:nx
        weight_sum = 0.0
        velocity_sum = 0.0
        for k in 1:nv
            w = Float64(cube[x, y, k])
            isfinite(w) || continue
            weight_sum += w
            velocity_sum += w * velocities_kms[k]
        end
        if weight_sum > 0
            out[x, y] = velocity_sum / weight_sum
        end
    end
    return out
end

column_density_map(cube, dv_kms::Real) =
    HI_COLUMN_DENSITY_FACTOR .* moment0(cube) .* Float64(dv_kms)

doppler_temperature_k(sigma_kms) = DOPPLER_TEMPERATURE_FACTOR .* sigma_kms .^ 2

function finite_weighted_mean(values, weights)
    total = 0.0
    weighted = 0.0
    for (v, w) in zip(values, weights)
        if isfinite(v) && isfinite(w) && w > 0
            total += w
            weighted += v * w
        end
    end
    return total > 0 ? weighted / total : NaN
end

function skewness(values)
    vals = [Float64(v) for v in values if isfinite(v)]
    isempty(vals) && return NaN
    m = mean(vals)
    s = std(vals)
    s == 0 && return NaN
    return mean(((vals .- m) ./ s).^3)
end

function rohsa_width_group(weighted_sigma_kms::Real)
    # Width bins are descriptive reporting groups for ROHSA products, not a
    # production VGT mask or a validated thermal phase classifier.
    return weighted_sigma_kms < 3 ? "narrow" :
           weighted_sigma_kms < 7 ? "intermediate" : "broad"
end

function rohsa_td_group(weighted_t_doppler_k::Real)
    # Doppler-temperature bins follow the usual HI phase shorthand and are used
    # only to write convenient grouped cubes and quick-look maps.
    return weighted_t_doppler_k < 500 ? "td_lt_500k" :
           weighted_t_doppler_k < 5000 ? "td_500_5000k" : "td_ge_5000k"
end

function rohsa_component_rows(region_label::AbstractString, n_gauss::Int, params,
                              components, velocities_kms::AbstractVector, dv_kms::Real)
    total = sum(c -> sum(Float64, c), components)
    rows = Dict{String, Any}[]
    for g in 1:n_gauss
        comp = components[g]
        amp = @view params[3g - 2, :, :]
        mu_chan = @view params[3g - 1, :, :]
        sig_chan = @view params[3g, :, :]
        moment = moment0(comp)
        sigma_kms = abs(Float64(dv_kms)) .* sig_chan
        # params hold zero-based centres (read_rohsa_parameters subtracts 1 from
        # ROHSA's 1-based mu), so mu = 0 is the first channel.
        mu_kms = begin
            first_velocity = first(velocities_kms)
            step = length(velocities_kms) > 1 ? velocities_kms[2] - velocities_kms[1] : Float64(dv_kms)
            first_velocity .+ step .* mu_chan
        end
        weighted_sigma = finite_weighted_mean(vec(sigma_kms), vec(moment))
        weighted_mu = finite_weighted_mean(vec(mu_kms), vec(moment))
        weighted_td = doppler_temperature_k(weighted_sigma)
        width_group = rohsa_width_group(weighted_sigma)
        td_group = rohsa_td_group(weighted_td)
        peak = maximum(amp)
        push!(rows, Dict{String, Any}(
            "region" => region_label,
            "n_gauss" => n_gauss,
            "component" => g,
            "width_group" => width_group,
            "td_group" => td_group,
            "inband_emission_fraction" => total > 0 ? sum(Float64, comp) / total : NaN,
            "weighted_sigma_kms" => weighted_sigma,
            "weighted_mu_kms" => weighted_mu,
            "weighted_t_doppler_k" => weighted_td,
            "amp_fill_fraction_gt_5pct_peak" => peak > 0 ? count(v -> v > 0.05peak, amp) / length(amp) : 0.0,
        ))
    end
    return rows
end

function rohsa_fit_metrics(region_label::AbstractString, n_gauss::Int, data_cube, model)
    residual = Float64.(model) .- Float64.(data_cube)
    data_moment = moment0(data_cube)
    model_moment = moment0(model)
    moment_residual = model_moment .- data_moment
    return Dict{String, Any}(
        "region" => region_label,
        "n_gauss" => n_gauss,
        "total_model_over_data" => sum(Float64, model) / sum(Float64, data_cube),
        "voxel_rmse" => sqrt(mean(residual .^ 2)),
        "voxel_mae" => mean(abs.(residual)),
        "voxel_residual_skew" => skewness(residual),
        "moment_relative_rmse" => sqrt(mean(moment_residual .^ 2)) / mean(data_moment),
        "moment_mean_relative_bias" => mean(moment_residual) / mean(data_moment),
        "moment_residual_skew" => skewness(moment_residual),
    )
end

function add_fraction_map!(dict, key::String, numerator, denominator)
    frac = fill(NaN, size(numerator))
    @inbounds for idx in eachindex(numerator, denominator)
        den = denominator[idx]
        if isfinite(den) && den > 0
            frac[idx] = numerator[idx] / den
        end
    end
    dict[key] = frac
    return dict
end

function write_rohsa_decomposition_products(crop_dir::AbstractString, crop_label::AbstractString,
                                            header, params, components, data_cube;
                                            n_gauss::Int=length(components),
                                            run_label::AbstractString="",
                                            fits_extension::AbstractString=".fits",
                                            write_cubes::Bool=true)
    nv = size(data_cube, 3)
    velocities_kms, dv_kms = velocity_axis_kms(header, nv)
    file_prefix = isempty(run_label) ? "rohsa_n$(n_gauss)" : run_label
    tree_label = isempty(run_label) ? "n$(n_gauss)" : run_label
    fits_ext = normalize_fits_extension(fits_extension)
    root = joinpath(crop_dir, "rohsa", tree_label)
    map_dir = joinpath(root, "maps")
    cube_dir = joinpath(root, "cubes")
    mkpath(map_dir)
    mkpath(cube_dir)

    model = zeros(Float32, size(components[1]))
    width_groups = Dict(
        "narrow" => zeros(Float32, size(components[1])),
        "intermediate" => zeros(Float32, size(components[1])),
        "broad" => zeros(Float32, size(components[1])),
    )
    temperature_groups = Dict(
        "td_lt_500k" => zeros(Float32, size(components[1])),
        "td_500_5000k" => zeros(Float32, size(components[1])),
        "td_ge_5000k" => zeros(Float32, size(components[1])),
    )

    component_rows = rohsa_component_rows(crop_label, n_gauss, params, components, velocities_kms, dv_kms)
    component_paths = Dict{Int, String}()
    cube_paths = Dict{String, String}()
    map_paths = Dict{String, String}()
    for g in 1:n_gauss
        comp = components[g]
        model .+= comp
        row = component_rows[g]
        # Each Gaussian contributes to exactly one width group and one
        # Doppler-temperature group based on its moment-weighted summaries.
        width_groups[string(row["width_group"])] .+= comp
        temperature_groups[string(row["td_group"])] .+= comp

        if write_cubes
            component_path = joinpath(crop_dir, "$(file_prefix)_component_g$(g)$(fits_ext)")
            write_fits(component_path, comp; header=header)
            component_paths[g] = component_path
            cube_paths["component_g$(g)"] = component_path
        else
            component_paths[g] = ""
        end

        amp = Float32.(@view params[3g - 2, :, :])
        mu_chan = @view params[3g - 1, :, :]
        sig_chan = @view params[3g, :, :]
        step = length(velocities_kms) > 1 ? velocities_kms[2] - velocities_kms[1] : Float64(dv_kms)
        mu_kms = Float32.(first(velocities_kms) .+ step .* mu_chan)  # zero-based centres
        sigma_kms = Float32.(abs(Float64(dv_kms)) .* sig_chan)
        t_doppler = Float32.(doppler_temperature_k.(Float64.(sigma_kms)))
        nhi = Float32.(column_density_map(comp, dv_kms))

        for (name, data) in (
            "amp_k" => amp,
            "mu_kms" => mu_kms,
            "sigma_kms" => sigma_kms,
            "t_doppler_k" => t_doppler,
            "nhi_cm2" => nhi,
        )
            path = joinpath(map_dir, "$(file_prefix)_component_g$(g)_$(name)$(fits_ext)")
            write_fits(path, data; header=header)
            map_paths["component_g$(g)_$(name)"] = path
        end
    end

    # Cubes are written once, at the region root, and only when asked for:
    # every one of them can be rebuilt from the ROHSA parameter file.
    model_path = write_cubes ? joinpath(crop_dir, "$(file_prefix)_model$(fits_ext)") : ""
    if write_cubes
        write_fits(model_path, model; header=header)
        cube_paths["model"] = model_path
    end

    group_paths = Dict{String, String}()
    group_nhi = Dict{String, Matrix{Float64}}()
    all_groups = merge(width_groups, temperature_groups)
    for (group, cube) in all_groups
        if write_cubes
            path = joinpath(crop_dir, "$(file_prefix)_$(group)_group$(fits_ext)")
            write_fits(path, cube; header=header)
            group_paths[group] = path
            cube_paths["$(group)_group"] = path
        end
        nhi = column_density_map(cube, dv_kms)
        group_nhi[group] = nhi
        map_path = joinpath(map_dir, "$(file_prefix)_$(group)_nhi_cm2$(fits_ext)")
        write_fits(map_path, Float32.(nhi); header=header)
        map_paths["$(group)_nhi_cm2"] = map_path
    end

    model_nhi = column_density_map(model, dv_kms)
    data_nhi = column_density_map(data_cube, dv_kms)
    residual = Float32.(Float64.(model) .- Float64.(data_cube))
    residual_moment = moment0(model) .- moment0(data_cube)
    residual_nhi = model_nhi .- data_nhi
    if write_cubes
        residual_cube_path = joinpath(cube_dir, "$(file_prefix)_residual_cube$(fits_ext)")
        write_fits(residual_cube_path, residual; header=header)
        cube_paths["residual"] = residual_cube_path
    end
    # Residual convention is model - data throughout the report and products.
    for (name, data) in (
        "model_nhi_cm2" => model_nhi,
        "data_nhi_cm2" => data_nhi,
        "residual_moment0_kkms" => residual_moment .* dv_kms,
        "residual_nhi_cm2" => residual_nhi,
    )
        path = joinpath(map_dir, "$(file_prefix)_$(name)$(fits_ext)")
        write_fits(path, Float32.(data); header=header)
        map_paths[name] = path
    end

    for group in ("narrow", "intermediate", "broad", "td_lt_500k", "td_500_5000k", "td_ge_5000k")
        # Fractions are normalized by the reconstructed ROHSA model column
        # density, so they describe the decomposition rather than raw data.
        add_fraction_map!(group_nhi, "$(group)_fraction", group_nhi[group], model_nhi)
        path = joinpath(map_dir, "$(file_prefix)_$(group)_fraction$(fits_ext)")
        write_fits(path, Float32.(group_nhi["$(group)_fraction"]); header=header)
        map_paths["$(group)_fraction"] = path
    end

    fit_row = rohsa_fit_metrics(crop_label, n_gauss, data_cube, model)
    return (;
        root,
        map_dir,
        cube_dir,
        velocities_kms,
        dv_kms,
        component_rows,
        fit_row,
        component_paths,
        cube_paths,
        map_paths,
        model_path,
        group_paths,
        model,
        groups=all_groups,
        model_nhi,
        data_nhi,
        group_nhi,
    )
end
