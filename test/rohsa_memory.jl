using Test

module RohsaMemoryTest
    include(joinpath(@__DIR__, "..", "rohsa", "memory.jl"))
end

@testset "rohsa_memory_estimate derivative term tracks f_g_cube_fast, not nv" begin
    m = RohsaMemoryTest

    small = m.rohsa_memory_estimate(16, 16, 100, 4, 10)
    large_nv = m.rohsa_memory_estimate(16, 16, 100_000, 4, 10)

    # The production gradient routine (f_g_cube_fast) accumulates into a
    # 3N x ny x nx map, independent of the number of spectral channels.
    # A regression that reintroduces the old f_g_cube nv-scaled tensor here
    # would make the launch guard wildly overestimate memory for large cubes.
    @test small.derivative_map == large_nv.derivative_map
    @test small.derivative_map == 8 * 3 * small.n_gauss * small.padded_side^2

    baseline = m.rohsa_memory_estimate(64, 64, 100, 6, 10)
    more_gauss = m.rohsa_memory_estimate(64, 64, 100, 12, 10)
    @test more_gauss.derivative_map == 2 * baseline.derivative_map
    @test more_gauss.peak > baseline.peak
end
