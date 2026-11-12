using Printf
using Statistics
using Test
using VGT
using RohsaLMC
using TOML

Base.include(@__MODULE__, joinpath(@__DIR__, "..", "rohsa", "configs.jl"))
Base.include(@__MODULE__, joinpath(@__DIR__, "..", "rohsa", "selection.jl"))
Base.include(@__MODULE__, joinpath(@__DIR__, "..", "rohsa", "auto_config.jl"))

function base_auto_opts(; enabled="true", max_candidates="24")
    return Dict{String, String}(
        "auto-config-enabled" => enabled,
        "auto-config-profile" => "physical_decomposition",
        "auto-config-max-candidates" => max_candidates,
        "n-gauss" => "2",
        "lambda-amp" => "100.0",
        "lambda-mu" => "100.0",
        "lambda-sig" => "100.0",
        "lambda-var-amp" => "0.0",
        "lambda-var-mu" => "0.0",
        "lambda-var-sig" => "0.0",
        "lb-sig" => "1.0",
        "ub-sig" => "10.0",
        "auto-config-small-component-fraction" => "0.05",
        "auto-config-broad-fraction-limit" => "0.60",
        "auto-config-sigma-bound-tolerance" => "0.05",
        "auto-config-duplicate-velocity-kms" => "2.0",
        "auto-config-duplicate-sigma-kms" => "1.0",
        "auto-config-physical-weight" => "1.0",
        "auto-config-stability-weight" => "0.05",
    )
end

@testset "ROHSA auto-config candidates" begin
    disabled = base_auto_opts(; enabled="false")
    candidates = resolved_rohsa_candidates(disabled)
    @test length(candidates) == 1
    @test candidates[1]["candidate-id"] == ""
    @test candidates[1]["n-gauss"] == "2"

    opts = base_auto_opts(; max_candidates="3")
    opts["auto-config-n-gauss"] = "2,4"
    opts["auto-config-lambda-amp"] = "10.0,100.0"
    opts["auto-config-lambda-mu"] = "10.0,100.0"
    candidates = resolved_rohsa_candidates(opts)
    @test length(candidates) == 3
    @test [candidate["candidate-id"] for candidate in candidates] == ["c001_n2", "c002_n2", "c003_n4"]
    @test candidates[1]["lambda-amp"] == "10.0"
    @test candidates[end]["lambda-mu"] == "100.0"
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
        opts = load_rohsa_config(cfg)
        @test opts["fits-extension"] == ".fits.fz"
    end
end

@testset "ROHSA physical-decomposition score prefers stable components" begin
    opts = base_auto_opts()
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
    stable = apply_auto_config_scores!(copy(row), stable_params, stable_products, opts)

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
    bad = apply_auto_config_scores!(copy(row), bad_params, bad_products, opts)

    @test stable["selection_score"] < bad["selection_score"]
    @test bad["physical_penalty"] > stable["physical_penalty"]
end
