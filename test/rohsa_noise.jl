using Test
using Random

module RohsaInputsTest
    using Astronomy
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
