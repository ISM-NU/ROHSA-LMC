using Test

module RohsaMemoryTest
    include(joinpath(@__DIR__, "..", "rohsa", "memory.jl"))
end

@testset "rohsa_memory_estimate: channels only enter through the data cube" begin
    m = RohsaMemoryTest

    small = m.rohsa_memory_estimate(16, 16, 100, 4, 10)
    large_nv = m.rohsa_memory_estimate(16, 16, 100_000, 4, 10)

    # f_g_cube_fast keeps no nv-sized array besides the data. A regression that
    # adds the old f_g_cube 3N x nv x ny x nx tensor (or a residual cube) here
    # would make the launch guard wildly overestimate memory for large cubes.
    data(e) = e.data_original + e.data_padded
    @test large_nv.lower_bound - small.lower_bound == data(large_nv) - data(small)

    baseline = m.rohsa_memory_estimate(64, 64, 100, 6, 10)
    more_gauss = m.rohsa_memory_estimate(64, 64, 100, 12, 10)
    @test more_gauss.grid_params == 2 * baseline.grid_params
    @test more_gauss.peak > baseline.peak
end
