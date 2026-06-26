using Printf

const GIB = 1024.0^3

function rohsa_padded_side(nx::Int, ny::Int)
    side = max(nx, ny)
    side <= 1 && return 1
    # Upstream ROHSA pads spatial grids to a square power-of-two hierarchy.
    return 2 ^ ceil(Int, log2(side))
end

function rohsa_memory_estimate(nx::Int, ny::Int, nv::Int, n_gauss::Int, m::Int)
    padded = rohsa_padded_side(nx, ny)
    original_pixels = nx * ny
    padded_pixels = padded * padded
    n_beta = 3 * n_gauss * padded_pixels + n_gauss

    # This is a conservative lower-bound model of the largest upstream Fortran
    # arrays, plus L-BFGS-B workspace. It is a launch guard, not a profiler.
    data_original = 8 * nv * original_pixels
    data_padded = 8 * nv * padded_pixels
    grid_params = 8 * 3 * n_gauss * padded_pixels
    bounds_and_beta = 8 * 5 * n_beta
    lbfgs_work = 8 * (2 * m * n_beta + 5 * n_beta + 11 * m * m + 8 * m)
    derivative_tensor = 8 * 3 * n_gauss * nv * padded_pixels
    residual = 8 * nv * padded_pixels
    small_maps = 8 * 20 * padded_pixels

    lower_bound = data_original + data_padded + grid_params + bounds_and_beta +
                  lbfgs_work + derivative_tensor + residual + small_maps
    peak = ceil(Int, lower_bound * 1.15)
    return (;
        nx,
        ny,
        nv,
        n_gauss,
        m,
        padded_side=padded,
        data_original,
        data_padded,
        grid_params,
        lbfgs_work,
        derivative_tensor,
        lower_bound,
        peak,
        peak_gib=peak / GIB,
    )
end

function memory_budget_gib(opts)
    fraction = parse(Float64, opts["safe-memory-fraction"])
    0 < fraction <= 1 ||
        throw(ArgumentError("--safe-memory-fraction must be in (0, 1], got $(opts["safe-memory-fraction"])"))
    raw = opts["memory-gb"]
    # Keep automatic budgets below physical RAM so the Fortran optimizer has
    # room for allocator overhead and the Julia wrapper.
    base = lowercase(string(raw)) == "auto" ? Sys.total_memory() / GIB : parse(Float64, raw)
    base > 0 || throw(ArgumentError("--memory-gb must be positive or auto"))
    return base * fraction
end

function suggested_spatial_stride(image_nx::Int, image_ny::Int, nv::Int, n_gauss::Int, m::Int, budget_gib::Real)
    # Return the smallest stride that fits the configured budget so exploratory
    # runs preserve as much spatial resolution as the preflight estimate allows.
    for stride in 1:128
        nx = cld(image_nx, stride)
        ny = cld(image_ny, stride)
        rohsa_memory_estimate(nx, ny, nv, n_gauss, m).peak_gib <= budget_gib && return stride
    end
    return 128
end

function format_gib(bytes::Integer)
    return @sprintf("%.1f GiB", bytes / GIB)
end
