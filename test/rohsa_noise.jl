using Test
using Random
using Statistics: median

module RohsaInputsTest
    using VGT
    using FITSIO
    using Printf
    using Statistics
    include(joinpath(@__DIR__, "..", "rohsa", "inputs.jl"))
end

@testset "ROHSA noise estimate ignores mask fill and clipping" begin
    rng = MersenneTwister(7)
    sigma_true = 0.3
    cube = Float32.(sigma_true .* randn(rng, 24, 24, 40))
    raw_sigma = RohsaInputsTest.robust_channel_difference_noise(cube)
    @test isapprox(raw_sigma, sigma_true; rtol=0.1)

    # Masked voxels (NaN and exact-zero fill) must not pull sigma toward zero.
    masked = copy(cube)
    masked[1:12, :, :] .= NaN32
    masked[13:20, :, :] .= 0.0f0
    @test isapprox(RohsaInputsTest.robust_channel_difference_noise(masked), sigma_true; rtol=0.15)
end

@testset "ROHSA spectral kernel: noise factor and line broadening" begin
    for name in ("hanning", "gauss1", "gauss2")
        w = RohsaInputsTest.spectral_kernel_weights(name)
        rng = MersenneTwister(11)
        noise = Float32.(randn(rng, 40, 40, 200))
        smoothed = RohsaInputsTest.smooth_spectra(noise, w)
        # Reported ROHSA noise = raw sigma × sqrt(Σk²) must match the smoothed cube's rms.
        @test isapprox(sqrt(sum(abs2, smoothed[:, :, 20:180]) / length(smoothed[:, :, 20:180])), RohsaInputsTest.spectral_kernel_noise_factor(w); rtol=0.03)

        # A Gaussian line widens by the kernel sigma in quadrature; area is preserved.
        v = collect(1.0:201.0)
        line = reshape(Float32.(exp.(-0.5 .* ((v .- 101) ./ 2.0) .^ 2)), 1, 1, :)
        out = vec(RohsaInputsTest.smooth_spectra(line, w))
        width = sqrt(sum(out .* (v .- 101) .^ 2) / sum(out))
        @test isapprox(width, hypot(2.0, RohsaInputsTest.spectral_kernel_sigma(w)); rtol=1e-3)
        @test isapprox(sum(out), sum(line); rtol=1e-4)
    end
    @test RohsaInputsTest.spectral_kernel_weights("none") == [1.0]
    @test_throws ArgumentError RohsaInputsTest.spectral_kernel_weights("box3")
    # NaN mask is preserved and does not leak into neighbours.
    c = ones(Float32, 1, 1, 9); c[1, 1, 5] = NaN32
    s = RohsaInputsTest.smooth_spectra(c, RohsaInputsTest.spectral_kernel_weights("hanning"))
    @test isnan(s[1, 1, 5]) && all(s[1, 1, [1:4; 6:9]] .≈ 1)
end

@testset "ROHSA input cleaning: NaN to 0, negatives kept unless clipped" begin
    c = RohsaInputsTest.clean_voxel
    @test c(NaN32, true) == 0.0f0 && c(NaN32, false) == 0.0f0
    @test c(Inf32, false) == 0.0f0
    @test c(-0.5f0, true) == 0.0f0          # historical default
    @test c(-0.5f0, false) == -0.5f0        # clip_negative = false keeps zero-mean noise
    @test c(1.25f0, true) == 1.25f0
    # Why it matters: clipping zero-mean noise leaves a positive pedestal.
    rng = MersenneTwister(3)
    n = Float32.(randn(rng, 100_000))
    @test sum(c.(n, true)) / length(n) > 0.35           # ~0.4 sigma
    @test abs(sum(c.(n, false)) / length(n)) < 0.02
end

@testset "ROHSA line-free noise map: per-pixel level, unbiased by lines and channel correlation" begin
    # Synthetic cube: noise level rising from 0.2 to 0.4 across x, channel-correlated ([b, 1, b], lag 1 = 0.4, as
    # after resampling), and one bright line per pixel (peak 30x the noise, sigma 3 channels) at a varying centre.
    rng = MersenneTwister(5)
    nx, ny, nv = 32, 24, 120
    level = [0.2 + 0.2 * (x - 1) / (nx - 1) for x in 1:nx, y in 1:ny]
    b = (1 - sqrt(1 - 2 * 0.4^2)) / (2 * 0.4)
    w = [b, 1.0, b] ./ sqrt(1 + 2b^2)
    white = randn(rng, nx, ny, nv + 2)
    cube = zeros(Float32, nx, ny, nv)
    for c in 1:nv, y in 1:ny, x in 1:nx
        mu = 40 + 40 * (x + y) / (nx + ny)
        cube[x, y, c] = level[x, y] * (w[1] * white[x, y, c] + w[2] * white[x, y, c + 1] + w[3] * white[x, y, c + 2]) +
                        30 * 0.3 * exp(-0.5 * ((c - mu) / 3)^2)
    end
    scalar = RohsaInputsTest.robust_channel_difference_noise(cube)
    m = RohsaInputsTest.line_free_noise_map(cube; smooth_px=5)
    # interior pixels (the 5 px box mixes levels at the x edges); per-channel truth = level
    ratio = m[4:end-3, :] ./ level[4:end-3, :]
    @test isapprox(median(ratio), 1.0; atol=0.05)
    @test all(0.85 .< ratio .< 1.15)
    @test m[end, 1] > 1.6 * m[1, 1]                       # it follows the level across the field
    @test abs(scalar / median(level) - 1) > 0.1           # the default estimator is off on this cube

    # Fixed line-free channels give the same level; ranges are 1-based and inclusive.
    ch = RohsaInputsTest.parse_channel_ranges("1-25,100-120", nv)
    @test ch == [1:25; 100:120]
    mc = RohsaInputsTest.line_free_noise_map(cube; channels=ch, smooth_px=5)
    @test isapprox(median(mc[4:end-3, :] ./ level[4:end-3, :]), 1.0; atol=0.06)
    @test RohsaInputsTest.parse_channel_ranges("", nv) === nothing
    @test_throws ArgumentError RohsaInputsTest.parse_channel_ranges("0-5", nv)

    # Mask fill (NaN, exact zero) is ignored; a pixel with no usable voxels is filled from its neighbours.
    c2 = copy(cube); c2[1:4, 1:4, :] .= NaN32; c2[10, 10, :] .= 0.0f0
    m2 = RohsaInputsTest.line_free_noise_map(c2; smooth_px=5)
    @test all(isfinite, m2) && all(m2 .> 0)
    @test isapprox(m2[10, 10], level[10, 10]; rtol=0.15)
end

@testset "ROHSA line-free noise map: emission over most of the band" begin
    # Wide bright emission (two broad lines covering about 70 % of the channels, as in Milky Way H I) must not set the
    # emission threshold.
    rng = MersenneTwister(9)
    nx, ny, nv = 16, 16, 106
    cube = Float32.(0.3 .* randn(rng, nx, ny, nv))
    for c in 1:nv
        cube[:, :, c] .+= Float32(15 * exp(-0.5 * ((c - 45) / 12)^2) + 8 * exp(-0.5 * ((c - 70) / 6)^2))
    end
    m = RohsaInputsTest.line_free_noise_map(cube; smooth_px=5)
    @test isapprox(median(m), 0.3; rtol=0.08)
    @test RohsaInputsTest.robust_channel_difference_noise(cube) > 0.33   # line slopes raise the default estimate
end

@testset "ROHSA live channel ranges and the range file" begin
    # Pixels of a mosaic: full band; channels 1-4 and 15-20 fill (exact zeros) with a negative live voxel at the first
    # live channel; all NaN (beyond the survey edge); all zeros.
    nx, ny, nv = 3, 2, 20
    cube = fill(0.5f0, nx, ny, nv)
    cube[2, 1, [1:4; 15:20]] .= 0.0f0
    cube[2, 1, 5] = -0.3f0
    cube[3, 1, :] .= NaN32
    cube[1, 2, :] .= 0.0f0
    first, last = RohsaInputsTest.live_channel_ranges(cube)
    @test (first[1, 1], last[1, 1]) == (1, 20)
    @test (first[2, 1], last[2, 1]) == (5, 14)
    @test first[3, 1] > last[3, 1] && first[1, 2] > last[1, 2]
    mktempdir() do dir
        path = RohsaInputsTest.write_range_map(joinpath(dir, "r.dat"), first, last)
        lines = readlines(path)
        @test lines[1] == "$ny\t$nx"                   # header "ny nx", rows "y x first last", 0-based
        @test lines[2] == "0\t0\t0\t19"
        @test lines[3] == "0\t1\t4\t13"
        @test lines[4] == "0\t2\t0\t-1"                # empty range
        @test lines[5] == "1\t0\t0\t-1"
        @test length(lines) == 1 + nx * ny
    end
end

@testset "ROHSA line-free noise map ignores the fill channels of a mosaic block" begin
    # Same noise everywhere; in the right half the first 40 channels are fill (exact zeros), as outside an observing
    # block's channel range. The map there must not change: zeros are never noise samples.
    rng = MersenneTwister(21)
    nx, ny, nv = 24, 16, 120
    cube = Float32.(0.25 .* randn(rng, nx, ny, nv))
    for c in 1:nv
        cube[:, :, c] .+= Float32(6 * exp(-0.5 * ((c - 70) / 4)^2))
    end
    zeroed = copy(cube)
    zeroed[13:end, :, 1:40] .= 0.0f0
    m = RohsaInputsTest.line_free_noise_map(cube; smooth_px=5)
    mz = RohsaInputsTest.line_free_noise_map(zeroed; smooth_px=5)
    @test isapprox(median(mz[15:end, :]), 0.25; rtol=0.08)
    @test all(0.85 .< mz[15:end, :] ./ m[15:end, :] .< 1.15)
end

@testset "ROHSA attenuation weight and reliable channel ranges (synthetic mosaic taper)" begin
    # Two observing blocks: A (x <= 98) covers channels 1-60, B (x >= 99) only 21-60 (channels 1-20 exact zeros). The
    # mosaic averaged B's emptiness into A's one-block channels 1-20, scaling signal and noise by
    # w(x) = (99 - x) / 60 (1 at x <= 39, 0.02 at x = 98). A bright line in channels 30-50 everywhere. With 16 px tiles
    # the tile centres have w = 1, 1, 0.98, 0.71, 0.44, 0.18 on the A side, and the seam crosses tile 7, which holds
    # only two columns of A (96 pixels, too few differences to measure: filled from the neighbours). Reliable range
    # (w >= 0.8): from channel 1 in the first three tiles, from channel 21 or 22 beyond (channel 21 forms a pair with
    # the attenuated channel 20 and can be dropped too: conservative by design); block B pixels keep their channel 21,
    # since the pair (20, 21) is outside their live range.
    rng = MersenneTwister(31)
    nx, ny, nv = 160, 48, 60
    taper(x) = clamp((99 - x) / 60, 0.0, 1.0)
    cube = zeros(Float32, nx, ny, nv)
    for c in 1:nv, y in 1:ny, x in 1:nx
        v = 0.2 * randn(rng) + 8 * exp(-0.5 * ((c - 40) / 3)^2) + 0.5 * exp(-0.5 * ((c - 12) / 4)^2)
        if c <= 20
            v = x >= 99 ? 0.0 : taper(x) * v
        end
        cube[x, y, c] = Float32(v)
    end
    first, last = RohsaInputsTest.live_channel_ranges(cube)
    @test all(first[1:98, :] .== 1) && all(first[99:end, :] .== 21) && all(last .== nv)
    w, ref = RohsaInputsTest.attenuation_weight(cube, first, last; tile=16, smooth=3)
    @test isapprox(ref, 0.2; rtol=0.05)
    # The weight follows the taper in the one-block channels (the tiles next to the seam have measured neighbours on
    # one side only, so the median filter pulls them towards the interior), and is >= 0.8 in the full channels: the
    # line raises the difference noise, it never lowers it.
    # w[:, :, k] is the pair of channels k and k + 1: pairs 1-19 attenuated, 21-59 full.
    @test size(w) == (10, 3, nv - 1)
    for tx in 1:5
        @test isapprox(median(w[tx, :, 1:19]), taper((tx - 1) * 16 + 8.5); atol=0.08)
    end
    @test median(w[6, :, 1:19]) < 0.5 && median(w[7, :, 1:19]) < 0.5
    @test all(w[:, :, 21:59] .>= 0.8)
    rf, rl = RohsaInputsTest.reliable_channel_ranges(first, last, w; tile=16, wmin=0.8)
    @test all(rl .== nv)
    @test all(rf[1:48, :] .== 1)
    @test all(21 .<= rf[49:98, :] .<= 22)        # attenuated channels dropped, also in the tile the seam crosses
    @test all(rf[99:end, :] .== 21)              # block B: its own live range
    # Pixels without data stay empty.
    c2 = copy(cube); c2[1:4, 1:4, :] .= NaN32
    f2, l2 = RohsaInputsTest.live_channel_ranges(c2)
    r2f, r2l = RohsaInputsTest.reliable_channel_ranges(f2, l2, w; tile=16)
    @test all(r2f[1:4, 1:4] .> r2l[1:4, 1:4])

    # The line-free noise map restricted to the reliable range is not biased low by the attenuated channels.
    m_all = RohsaInputsTest.line_free_noise_map(cube; smooth_px=5)
    m_rel = RohsaInputsTest.line_free_noise_map(cube; smooth_px=5, ranges=(rf, rl))
    @test median(m_all[70:92, :]) < 0.17
    @test isapprox(median(m_rel[70:92, :]), 0.2; rtol=0.1)
end
