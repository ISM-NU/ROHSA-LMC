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
