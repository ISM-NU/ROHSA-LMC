using Astronomy
using Test

@testset "ROHSA velocity layers" begin
    for offset in (-60.0, -20.0, -5.0, 0.0, 12.0, 20.0, 45.0)
        w = Astronomy.velocity_layer_weights(offset)
        @test sum(w) ≈ 1
        @test all(>=(0), w)
    end
    @test Astronomy.velocity_layer_weights(-40.0)[1] > 0.99
    @test Astronomy.velocity_layer_weights(0.0)[2] > 0.99
    @test Astronomy.velocity_layer_weights(40.0)[3] > 0.99

    # Two pixels, two Gaussians: one near the reference, one 40 km/s below it.
    nx, ny, nv = 2, 1, 120
    params = zeros(3 * 2, nx, ny)
    for x in 1:nx
        params[1:3, x, 1] .= (10.0, 60.0, 4.0)
        params[4:6, x, 1] .= (5.0, 20.0, 6.0)
    end
    reference = fill(60.0, nx, ny)
    cubes = Astronomy.build_layer_cubes(params, nv, 1.0; reference)
    @test cubes[:low] .+ cubes[:main] .+ cubes[:high] ≈ cubes[:model]
    @test sum(cubes[:low]) / sum(cubes[:model]) ≈ 5 * 6 / (10 * 4 + 5 * 6) rtol = 1e-3
    @test maximum(cubes[:high]) < 2e-3 * maximum(cubes[:model])  # soft-edge leakage only

    ref = Astronomy.reference_centroid_channels(params; sigma_pixels=1.0)
    @test all(isfinite, ref)
    @test 20 < ref[1, 1] < 60

    # Phases by width: σ = 1 km/s (121 K) cold, 4 km/s (1936 K) unstable, 8 km/s warm.
    @test Astronomy.gaussian_phase(1.0, 5.0) === :cold
    @test Astronomy.gaussian_phase(1.0, 0.5; cold_min_amp=1.0) === :none
    @test Astronomy.gaussian_phase(4.0, 5.0) === :unstable
    @test Astronomy.gaussian_phase(8.0, 5.0) === :warm
    pp = zeros(9, 1, 1)
    pp[:, 1, 1] .= (10.0, 40.0, 1.0, 6.0, 50.0, 4.0, 3.0, 70.0, 8.0)
    ph = Astronomy.build_phase_cubes(pp, 120, 1.0)
    @test ph[:cold] .+ ph[:unstable] .+ ph[:warm] ≈ ph[:model]
    @test sum(ph[:cold]) / sum(ph[:model]) ≈ 10 * 1 / (10 * 1 + 6 * 4 + 3 * 8) rtol = 1e-3
    ph2 = Astronomy.build_phase_cubes(pp, 120, 1.0; cold_min_amp=20.0)
    @test sum(ph2[:cold]) == 0
end
