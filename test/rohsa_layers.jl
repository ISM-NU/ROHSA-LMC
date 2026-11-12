using Astronomy
using RohsaLMC
using Test

@testset "ROHSA velocity layers" begin
    for offset in (-60.0, -20.0, -5.0, 0.0, 12.0, 20.0, 45.0)
        w = RohsaLMC.velocity_layer_weights(offset)
        @test sum(w) ≈ 1
        @test all(>=(0), w)
    end
    @test RohsaLMC.velocity_layer_weights(-40.0)[1] > 0.99
    @test RohsaLMC.velocity_layer_weights(0.0)[2] > 0.99
    @test RohsaLMC.velocity_layer_weights(40.0)[3] > 0.99

    # Two pixels, two Gaussians: one near the reference, one 40 km/s below it.
    nx, ny, nv = 2, 1, 120
    params = zeros(3 * 2, nx, ny)
    for x in 1:nx
        params[1:3, x, 1] .= (10.0, 60.0, 4.0)
        params[4:6, x, 1] .= (5.0, 20.0, 6.0)
    end
    reference = fill(60.0, nx, ny)
    cubes = RohsaLMC.build_layer_cubes(params, nv, 1.0; reference)
    @test cubes[:low] .+ cubes[:main] .+ cubes[:high] ≈ cubes[:model]
    @test sum(cubes[:low]) / sum(cubes[:model]) ≈ 5 * 6 / (10 * 4 + 5 * 6) rtol = 1e-3
    @test maximum(cubes[:high]) < 2e-3 * maximum(cubes[:model])  # soft-edge leakage only

    ref = RohsaLMC.reference_centroid_channels(params; sigma_pixels=1.0)
    @test all(isfinite, ref)
    @test 20 < ref[1, 1] < 60

    # Phases by width: σ = 1 km/s (121 K) cold, 4 km/s (1936 K) unstable, 12 km/s (17000 K) warm.
    for (sig, k) in ((1.0, 1), (4.0, 2), (12.0, 3))
        w = RohsaLMC.gaussian_phase_weights(sig, 5.0)
        @test sum(w) ≈ 1
        @test w[k] > 0.99
    end
    @test RohsaLMC.gaussian_phase_weights(1.0, 0.2; cold_min_amp=1.0)[1] < 0.05
    @test 0.4 < RohsaLMC.gaussian_phase_weights(sqrt(500 / 121.1), 5.0)[1] < 0.6  # soft edge at 500 K
    pp = zeros(9, 1, 1)
    pp[:, 1, 1] .= (10.0, 40.0, 1.0, 6.0, 50.0, 4.0, 3.0, 70.0, 8.0)
    ph = RohsaLMC.build_phase_cubes(pp, 120, 1.0)
    @test ph[:cold] .+ ph[:unstable] .+ ph[:warm] ≈ ph[:model] rtol = 1e-4
    @test sum(ph[:cold]) / sum(ph[:model]) ≈ 10 * 1 / (10 * 1 + 6 * 4 + 3 * 8) rtol = 1e-2
    ph2 = RohsaLMC.build_phase_cubes(pp, 120, 1.0; cold_min_amp=100.0)
    @test sum(ph2[:cold]) < 1e-2 * sum(ph2[:model])
    # Edges given as FWHM: a Gaussian exactly at the edge is half narrow; equal edges leave no unstable gas.
    t7 = RohsaLMC.fwhm_doppler_temperature(7.0)
    @test 0.4 < RohsaLMC.gaussian_phase_weights(7.0 / 2.3548, 5.0; cold_max_k=t7, warm_min_k=t7)[1] < 0.6
    ph3 = RohsaLMC.build_phase_cubes(pp, 120, 1.0; cold_max_k=t7, warm_min_k=t7)
    @test sum(ph3[:unstable]) < 1e-6 * sum(ph3[:model])
    @test ph3[:cold] .+ ph3[:warm] ≈ ph3[:model] rtol = 1e-4
end
