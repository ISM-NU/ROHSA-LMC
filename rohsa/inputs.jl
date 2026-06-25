const MAD_NORMAL = 0.6744897501960817

Base.@kwdef struct InputRegion
    label::String
    x0::Int
    y0::Int
    nx::Int
    ny::Int
    source_nx::Int = nx
    source_ny::Int = ny
    x_stride::Int = 1
    y_stride::Int = 1
end

function region_from_options(input::AbstractString, image_nx::Int, image_ny::Int, opts)
    crop_size = parse(Int, opts["crop-size"])
    spatial_stride = parse(Int, opts["spatial-stride"])
    spatial_stride >= 1 ||
        throw(ArgumentError("--spatial-stride must be >= 1, got $spatial_stride"))
    if crop_size <= 0
        nx = cld(image_nx, spatial_stride)
        ny = cld(image_ny, spatial_stride)
        return InputRegion(
            label=spatial_stride == 1 ? "full" : "full_stride$(spatial_stride)",
            x0=1,
            y0=1,
            nx=nx,
            ny=ny,
            source_nx=image_nx,
            source_ny=image_ny,
            x_stride=spatial_stride,
            y_stride=spatial_stride,
        )
    end

    crop_size <= min(image_nx, image_ny) ||
        throw(ArgumentError("--crop-size must fit inside input image dimensions ($image_nx, $image_ny)"))
    x0 = parse(Int, opts["crop-x0"])
    y0 = parse(Int, opts["crop-y0"])
    1 <= x0 <= image_nx - crop_size + 1 ||
        throw(ArgumentError("--crop-x0=$x0 does not fit crop size $crop_size inside image width $image_nx"))
    1 <= y0 <= image_ny - crop_size + 1 ||
        throw(ArgumentError("--crop-y0=$y0 does not fit crop size $crop_size inside image height $image_ny"))

    return InputRegion(
        label=spatial_stride == 1 ? @sprintf("crop_y%d_x%d", y0, x0) :
              @sprintf("crop_y%d_x%d_stride%d", y0, x0, spatial_stride),
        x0=x0,
        y0=y0,
        nx=cld(crop_size, spatial_stride),
        ny=cld(crop_size, spatial_stride),
        source_nx=crop_size,
        source_ny=crop_size,
        x_stride=spatial_stride,
        y_stride=spatial_stride,
    )
end

function maybe_reused_region(output_dir::AbstractString, opts, image_nx::Int, image_ny::Int)
    manifest_path = joinpath(output_dir, "manifest.toml")
    region_path = joinpath(output_dir, "region_selection.csv")
    (isfile(manifest_path) && isfile(region_path)) || return nothing
    manifest = TOML.parsefile(manifest_path)
    old_opts = get(manifest, "options", Dict{String, Any}())
    for key in ("input", "crop-size", "crop-x0", "crop-y0", "spatial-stride")
        string(get(old_opts, key, "")) == string(opts[key]) || return nothing
    end
    lines = readlines(region_path)
    length(lines) >= 2 || return nothing
    header = split(lines[1], ",")
    col(name) = something(findfirst(==(name), header), 0)
    needed = ["region", "x0", "y0", "nx", "ny"]
    all(col(name) > 0 for name in needed) || return nothing
    parts = split(lines[2], ",")
    region = InputRegion(
        label=parts[col("region")],
        x0=parse(Int, parts[col("x0")]),
        y0=parse(Int, parts[col("y0")]),
        nx=parse(Int, parts[col("nx")]),
        ny=parse(Int, parts[col("ny")]),
        source_nx=parse(Int, parts[col("nx")]),
        source_ny=parse(Int, parts[col("ny")]),
        x_stride=parse(Int, opts["spatial-stride"]),
        y_stride=parse(Int, opts["spatial-stride"]),
    )
    region.x0 + (region.nx - 1) * region.x_stride <= image_nx || return nothing
    region.y0 + (region.ny - 1) * region.y_stride <= image_ny || return nothing
    return region
end

function select_input_regions(input::AbstractString, image_nx::Int, image_ny::Int, opts; reuse::Bool)
    reused = reuse ? maybe_reused_region(opts["output-dir"], opts, image_nx, image_ny) : nothing
    reused !== nothing && return [reused]
    return [region_from_options(input, image_nx, image_ny, opts)]
end

function update_region_header!(header, region::InputRegion)
    header["CRPIX1"] = Float64(header["CRPIX1"]) - (region.x0 - 1)
    header["CRPIX2"] = Float64(header["CRPIX2"]) - (region.y0 - 1)
    if region.x_stride != 1
        header["CRPIX1"] = 1.0 + (Float64(header["CRPIX1"]) - 1.0) / region.x_stride
        header["CDELT1"] = Float64(header["CDELT1"]) * region.x_stride
    end
    if region.y_stride != 1
        header["CRPIX2"] = 1.0 + (Float64(header["CRPIX2"]) - 1.0) / region.y_stride
        header["CDELT2"] = Float64(header["CDELT2"]) * region.y_stride
    end
    return header
end

function robust_channel_difference_noise(cube)
    channel_diffs = vec(diff(Float64.(cube); dims=3))
    vals = [v for v in channel_diffs if isfinite(v)]
    isempty(vals) && return NaN
    med = median(vals)
    sigma = median(abs.(vals .- med)) / MAD_NORMAL / sqrt(2)
    return sigma > 0 ? sigma : NaN
end

function read_region_cube(input::AbstractString, region::InputRegion)
    f = FITS(input)
    try
        _, _, nv = size(f[1])
        x_range = region.x0:region.x_stride:region.x0 + (region.nx - 1) * region.x_stride
        y_range = region.y0:region.y_stride:region.y0 + (region.ny - 1) * region.y_stride
        return Float32.(read(f[1], x_range, y_range, 1:nv))
    finally
        close(f)
    end
end

function write_cube_dat(dat_path::AbstractString, cube)
    nx, ny, nv = size(cube)
    open(dat_path, "w") do io
        println(io, "$(nv)\t$(ny)\t$(nx)")
        for y in 1:ny, x in 1:nx, k in 1:nv
            @printf(io, "%d\t%d\t%d\t%.16f\n", k - 1, y - 1, x - 1, cube[x, y, k])
        end
    end
    return dat_path
end

function estimate_streamed_noise(input::AbstractString, region::InputRegion; max_samples::Int=2_000_000)
    f = FITS(input)
    vals = Float64[]
    try
        _, _, nv = size(f[1])
        total_diffs = region.nx * region.ny * max(nv - 1, 0)
        sample_step = max(1, cld(total_diffs, max_samples))
        counter = 0
        sizehint!(vals, min(total_diffs, max_samples))
        for y_out in 1:region.ny
            y_src = region.y0 + (y_out - 1) * region.y_stride
            row = Float32.(read(f[1], region.x0:region.x_stride:region.x0 + (region.nx - 1) * region.x_stride,
                                y_src:y_src, 1:nv))
            for x in 1:region.nx, k in 1:(nv - 1)
                counter += 1
                counter % sample_step == 0 || continue
                a = row[x, 1, k]
                b = row[x, 1, k + 1]
                if isfinite(a) && isfinite(b)
                    push!(vals, Float64(b) - Float64(a))
                end
            end
        end
    finally
        close(f)
    end
    isempty(vals) && return NaN
    med = median(vals)
    sigma = median(abs.(vals .- med)) / MAD_NORMAL / sqrt(2)
    return sigma > 0 ? sigma : NaN
end

function write_streamed_cube_dat(input::AbstractString, dat_path::AbstractString, region::InputRegion)
    f = FITS(input)
    try
        _, _, nv = size(f[1])
        open(dat_path, "w") do io
            println(io, "$(nv)\t$(region.ny)\t$(region.nx)")
            x_range = region.x0:region.x_stride:region.x0 + (region.nx - 1) * region.x_stride
            for y_out in 1:region.ny
                y_src = region.y0 + (y_out - 1) * region.y_stride
                row = Float32.(read(f[1], x_range, y_src:y_src, 1:nv))
                for x in 1:region.nx, k in 1:nv
                    v = row[x, 1, k]
                    clean = isfinite(v) && v > 0 ? v : 0.0f0
                    @printf(io, "%d\t%d\t%d\t%.16f\n", k - 1, y_out - 1, x - 1, clean)
                end
            end
        end
    finally
        close(f)
    end
    return dat_path
end

function write_rms_map(rms_path::AbstractString, nx::Int, ny::Int, sigma::Real)
    open(rms_path, "w") do io
        println(io, "$(ny)\t$(nx)")
        for y in 1:ny, x in 1:nx
            @printf(io, "%d\t%d\t%.16f\n", y - 1, x - 1, sigma)
        end
    end
    return rms_path
end

function write_region_products(input::AbstractString, header, region::InputRegion, region_dir::AbstractString;
                               retain_cube::Bool=true, write_region_fits::Bool=true)
    mkpath(region_dir)
    region_header = deepcopy(header)
    update_region_header!(region_header, region)
    fits_path = joinpath(region_dir, string(region.label, ".fits"))
    dat_path = joinpath(region_dir, string(region.label, ".dat"))

    cube = retain_cube || write_region_fits ? read_region_cube(input, region) : nothing
    if cube !== nothing
        cube .= map(v -> isfinite(v) && v > 0 ? v : 0.0f0, cube)
        write_region_fits && Astronomy.write_fits(fits_path, cube; header=region_header)
        write_cube_dat(dat_path, cube)
        sigma = robust_channel_difference_noise(cube)
    else
        write_streamed_cube_dat(input, dat_path, region)
        sigma = estimate_streamed_noise(input, region)
    end
    isfinite(sigma) || (sigma = 1.0)
    rms_path = joinpath(region_dir, "rms_map.dat")
    write_rms_map(rms_path, region.nx, region.ny, sigma)
    return (; fits_path, dat_path, rms_path, sigma, header=region_header, cube)
end
