using Printf
using Statistics
using Test
using VGT
using RohsaLMC
using TOML

Base.include(@__MODULE__, joinpath(@__DIR__, "..", "rohsa", "configs.jl"))
Base.include(@__MODULE__, joinpath(@__DIR__, "..", "rohsa", "selection.jl"))
Base.include(@__MODULE__, joinpath(@__DIR__, "..", "rohsa", "auto_config.jl"))

# A complete config file (every required key) loaded through load_rohsa_config.
const BASE_CONFIG_TOML = """
[input]
input = "cube.fits"
output_dir = "out"
crop_size = 0
crop_x0 = 1
crop_y0 = 1
reuse = true
write_visuals = false

[rohsa]
n_gauss = [2]
lambda_amp = 100.0
lambda_mu = 100.0
lambda_sig = 100.0
lambda_var_amp = 0.0
lambda_var_mu = 0.0
lambda_var_sig = 0.0
amp_fact_init = 0.66
sig_init = 4.0
lb_sig_init = 1.0
ub_sig_init = 12.0
lb_sig = 1.0
ub_sig = 10.0
init_option = "mean"
maxiter_init = 15000
maxiter = 800
m = 10
noise = false
regul = true
descent = true
lstd = 1
ustd = 20
iprint = -1
iprint_init = -1
save_grid = false
"""

function load_test_config(extra=""; enabled="true", max_candidates=24)
    mktempdir() do dir
        path = joinpath(dir, "rohsa.toml")
        write(path, BASE_CONFIG_TOML * """

        [auto_config]
        enabled = $enabled
        max_candidates = $max_candidates
        small_component_fraction = 0.05
        """ * extra)
        return load_rohsa_config(path)
    end
end

@testset "ROHSA auto-config candidates" begin
    disabled = resolved_rohsa_candidates(load_test_config(; enabled="false"))
    @test length(disabled) == 1
    @test disabled[1].id == ""
    @test disabled[1].n_gauss == 2

    cfg = load_test_config("""
    n_gauss = [2, 4]
    lambda_amp = [10.0, 100.0]
    lambda_mu = [10.0, 100.0]
    """; max_candidates=3)
    candidates = resolved_rohsa_candidates(cfg)
    @test length(candidates) == 3
    @test [candidate.id for candidate in candidates] == ["c001_n2", "c002_n2", "c003_n4"]
    @test candidates[1].lambda_amp == 10.0
    @test candidates[end].lambda_mu == 100.0

    # An axis the config leaves out holds the single [rohsa] value.
    @test all(c.lambda_sig == 100.0 for c in candidates)

    empty_axis = load_test_config("n_gauss = []\n")
    @test_throws ArgumentError resolved_rohsa_candidates(empty_axis)
end

@testset "ROHSA config normalizes FITS extension" begin
    mktempdir() do dir
        cfg = joinpath(dir, "rohsa.toml")
        open(cfg, "w") do io
            println(io, """
            [input]
            input = "cube.fits"
            output_dir = "out"
            crop_size = 0
            crop_x0 = 1
            crop_y0 = 1
            spatial_stride = 1
            reuse = true
            write_visuals = false
            derive_products = true
            write_region_fits = true
            fits_extension = ".fz"
            preflight_only = false
            force_full_cube = false
            memory_gb = "auto"
            safe_memory_fraction = 0.70

            [rohsa]
            n_gauss = [2]
            lambda_amp = 100.0
            lambda_mu = 100.0
            lambda_sig = 100.0
            lambda_var_amp = 0.0
            lambda_var_mu = 0.0
            lambda_var_sig = 0.0
            amp_fact_init = 0.66
            sig_init = 4.0
            lb_sig_init = 1.0
            ub_sig_init = 12.0
            lb_sig = 1.0
            ub_sig = 100.0
            init_option = "mean"
            maxiter_init = 15000
            maxiter = 800
            m = 10
            noise = false
            regul = true
            descent = true
            lstd = 1
            ustd = 20
            iprint = -1
            iprint_init = -1
            save_grid = false
            """)
        end
        @test load_rohsa_config(cfg).fits_extension == ".fits.fz"
    end
end

@testset "ROHSA physical-decomposition score prefers stable components" begin
    cfg = load_test_config()
    row = Dict{String, Any}(
        "total_model_over_data" => 1.0,
        "moment_relative_rmse" => 0.05,
        "moment_mean_relative_bias" => 0.0,
        "moment_residual_skew" => 0.0,
    )

    stable_params = zeros(Float64, 6, 4, 4)
    stable_params[1, :, :] .= 3.0
    stable_params[2, :, :] .= 2.0
    stable_params[3, :, :] .= 3.0
    stable_params[4, :, :] .= 2.0
    stable_params[5, :, :] .= 7.0
    stable_params[6, :, :] .= 4.0
    stable_products = (;
        component_rows=[
            Dict{String, Any}("inband_emission_fraction" => 0.60, "weighted_mu_kms" => 2.0, "weighted_sigma_kms" => 3.0),
            Dict{String, Any}("inband_emission_fraction" => 0.40, "weighted_mu_kms" => 7.0, "weighted_sigma_kms" => 4.0),
        ],
        group_nhi=Dict("broad_fraction" => fill(0.20, 4, 4)),
    )
    stable = apply_auto_config_scores!(copy(row), stable_params, stable_products, cfg)

    bad_params = copy(stable_params)
    bad_params[3, :, :] .= 1.0
    bad_params[6, :, :] .= 10.0
    bad_params[1, 2, 2] = 30.0
    bad_products = (;
        component_rows=[
            Dict{String, Any}("inband_emission_fraction" => 0.98, "weighted_mu_kms" => 2.0, "weighted_sigma_kms" => 3.0),
            Dict{String, Any}("inband_emission_fraction" => 0.02, "weighted_mu_kms" => 2.5, "weighted_sigma_kms" => 3.2),
        ],
        group_nhi=Dict("broad_fraction" => fill(0.90, 4, 4)),
    )
    bad = apply_auto_config_scores!(copy(row), bad_params, bad_products, cfg)

    @test stable["selection_score"] < bad["selection_score"]
    @test bad["physical_penalty"] > stable["physical_penalty"]
end
