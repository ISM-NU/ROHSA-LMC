function finite_values(arr)
    return [Float64(v) for v in arr if isfinite(v)]
end

function safe_colorrange(arr; symmetric=false, default=(0.0, 1.0))
    vals = finite_values(arr)
    isempty(vals) && return default
    if symmetric
        lim = maximum(abs, vals)
        lim = lim > 0 ? lim : 1.0
        return (-lim, lim)
    end
    lo = quantile(vals, 0.02)
    hi = quantile(vals, 0.98)
    if !(isfinite(lo) && isfinite(hi)) || lo == hi
        lo = minimum(vals)
        hi = maximum(vals)
    end
    if lo == hi
        pad = lo == 0 ? 1.0 : abs(lo) * 0.05
        return (lo - pad, hi + pad)
    end
    return (lo, hi)
end

moment0(cube) = dropdims(sum(Float64.(cube); dims=3); dims=3)

function heatmap_panel!(fig, row::Int, col::Int, data, title;
                        colormap=:viridis, colorrange=nothing)
    ax = CairoMakie.Axis(fig[row, col]; title=title, aspect=CairoMakie.DataAspect())
    CairoMakie.hidedecorations!(ax; grid=false)
    hm = CairoMakie.heatmap!(ax, data; colormap=colormap,
                             colorrange=colorrange === nothing ? safe_colorrange(data) : colorrange)
    CairoMakie.Colorbar(fig[row, col + 1], hm)
    return ax
end

function write_rohsa_visuals(region_dir, n_gauss::Int, data_cube, products; run_label::AbstractString="")
    visual_dir = joinpath(region_dir, "visuals")
    mkpath(visual_dir)
    file_prefix = isempty(run_label) ? "rohsa_n$(n_gauss)" : run_label
    data_moment = moment0(data_cube)
    model_moment = moment0(products.model)
    residual = model_moment .- data_moment
    group_moments = Dict(group => moment0(products.groups[group]) for group in keys(products.groups))
    shared_range = safe_colorrange(vcat(vec(data_moment), vec(model_moment)))
    width_groups = ("narrow", "intermediate", "broad")
    temp_groups = ("td_lt_500k", "td_500_5000k", "td_ge_5000k")
    group_range = safe_colorrange(vcat([vec(group_moments[group]) for group in width_groups]...))
    temp_range = safe_colorrange(vcat([vec(products.group_nhi[group]) for group in temp_groups]...))
    width_nhi_range = safe_colorrange(vcat([vec(products.group_nhi[group]) for group in width_groups]...))

    fig = CairoMakie.Figure(size=(1500, 900))
    heatmap_panel!(fig, 1, 1, data_moment, "data moment0"; colormap=:viridis,
                   colorrange=shared_range)
    heatmap_panel!(fig, 1, 3, model_moment, "ROHSA N=$(n_gauss) model moment0"; colormap=:viridis,
                   colorrange=shared_range)
    heatmap_panel!(fig, 1, 5, residual, "model - data moment0"; colormap=:balance,
                   colorrange=safe_colorrange(residual; symmetric=true))
    for (col, group) in zip((1, 3, 5), ("narrow", "intermediate", "broad"))
        heatmap_panel!(fig, 2, col, group_moments[group], "$(group) group moment0"; colormap=:magma,
                       colorrange=group_range)
    end
    path = joinpath(visual_dir, "$(file_prefix)_moment_panels.png")
    CairoMakie.save(path, fig)

    phase_fig = CairoMakie.Figure(size=(1500, 900))
    for (col, group) in zip((1, 3, 5), width_groups)
        heatmap_panel!(phase_fig, 1, col, products.group_nhi[group], "$(group) NHI"; colormap=:magma,
                       colorrange=width_nhi_range)
        heatmap_panel!(phase_fig, 2, col, products.group_nhi["$(group)_fraction"], "$(group) fraction";
                       colormap=:viridis, colorrange=(0.0, 1.0))
    end
    phase_path = joinpath(visual_dir, "$(file_prefix)_width_phase_maps.png")
    CairoMakie.save(phase_path, phase_fig)

    temp_fig = CairoMakie.Figure(size=(1500, 900))
    for (col, group) in zip((1, 3, 5), temp_groups)
        heatmap_panel!(temp_fig, 1, col, products.group_nhi[group], "$(group) NHI"; colormap=:magma,
                       colorrange=temp_range)
        heatmap_panel!(temp_fig, 2, col, products.group_nhi["$(group)_fraction"], "$(group) fraction";
                       colormap=:viridis, colorrange=(0.0, 1.0))
    end
    temp_path = joinpath(visual_dir, "$(file_prefix)_temperature_phase_maps.png")
    CairoMakie.save(temp_path, temp_fig)

    sv_fig = CairoMakie.Figure(size=(900, 620))
    ax = CairoMakie.Axis(sv_fig[1, 1];
                         title="ROHSA N=$(n_gauss) sigma-velocity diagram",
                         xlabel="weighted centroid velocity (km/s)",
                         ylabel="weighted sigma (km/s)")
    rows = sort(products.component_rows; by=row -> Int(row["component"]))
    xs = Float64[row["weighted_mu_kms"] for row in rows]
    ys = Float64[row["weighted_sigma_kms"] for row in rows]
    sizes = [max(8.0, 60.0 * sqrt(max(Float64(row["inband_emission_fraction"]), 0.0))) for row in rows]
    colors = [string(row["td_group"]) == "td_lt_500k" ? :dodgerblue3 :
              string(row["td_group"]) == "td_500_5000k" ? :goldenrod2 : :indianred3 for row in rows]
    CairoMakie.scatter!(ax, xs, ys; markersize=sizes, color=colors, strokecolor=:black, strokewidth=0.8)
    for row in rows
        CairoMakie.text!(ax, Float64(row["weighted_mu_kms"]), Float64(row["weighted_sigma_kms"]);
                         text=string("G", row["component"]), fontsize=12, align=(:left, :bottom))
    end
    CairoMakie.hlines!(ax, [sqrt(500 / Astronomy.DOPPLER_TEMPERATURE_FACTOR),
                            sqrt(5000 / Astronomy.DOPPLER_TEMPERATURE_FACTOR)];
                       color=[:dodgerblue3, :indianred3], linestyle=:dash)
    sv_path = joinpath(visual_dir, "$(file_prefix)_sigma_velocity.png")
    CairoMakie.save(sv_path, sv_fig)

    residual_fig = CairoMakie.Figure(size=(900, 520))
    residual_vals = finite_values(residual)
    axr = CairoMakie.Axis(residual_fig[1, 1];
                          title="ROHSA N=$(n_gauss) residual moment0 PDF",
                          xlabel="model - data moment0",
                          ylabel="count")
    !isempty(residual_vals) && CairoMakie.hist!(axr, residual_vals; bins=60, color=:gray65, strokecolor=:gray25)
    CairoMakie.vlines!(axr, [0.0]; color=:black, linestyle=:dash)
    residual_path = joinpath(visual_dir, "$(file_prefix)_residual_pdf.png")
    CairoMakie.save(residual_path, residual_fig)

    return [path, phase_path, temp_path, sv_path, residual_path]
end

function html_escape(text)
    return replace(string(text), "&" => "&amp;", "<" => "&lt;", ">" => "&gt;", "\"" => "&quot;")
end

function write_visual_index(path::AbstractString, visual_paths, output_paths)
    mkpath(dirname(path))
    base = dirname(path)
    open(path, "w") do io
        println(io, "<!doctype html><meta charset=\"utf-8\"><title>ROHSA preprocessing visuals</title>")
        println(io, "<style>body{font:14px system-ui,sans-serif;margin:24px;max-width:1200px}img{max-width:100%;border:1px solid #ccc;margin:8px 0 28px}code{background:#eee;padding:2px 4px}</style>")
        println(io, "<h1>ROHSA Preprocessing Visuals</h1>")
        println(io, "<h2>Figures</h2>")
        for path in visual_paths
            rel = relpath(path, base)
            println(io, "<figure><a href=\"", html_escape(rel), "\"><img src=\"", html_escape(rel), "\"></a><figcaption><code>", html_escape(rel), "</code></figcaption></figure>")
        end
        println(io, "<h2>Tables and Report</h2><ul>")
        for path in output_paths
            rel = relpath(path, base)
            println(io, "<li><a href=\"", html_escape(rel), "\"><code>", html_escape(rel), "</code></a></li>")
        end
        println(io, "</ul>")
    end
    return path
end
