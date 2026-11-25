using FITSIO
using Test
using Random

@testset "ROHSA preprocessing wrapper parses focused options" begin
    wrapper = Module(:RohsaComponentPipelineSmoke)
    Core.eval(wrapper, :(using Test))
    Base.include(wrapper, joinpath(@__DIR__, "..", "rohsa", "cli", "preprocess.jl"))

    cfg = wrapper.parse_args(String[])
    @test cfg.input == "data/exp_raw/Tb_combined_mw_full_hnv_resOG_reprojected_on_gass.fits"
    @test cfg.output_dir == "validation_output/rohsa_preprocessing"
    @test cfg.n_gauss == [6]
    @test cfg.write_visuals == true
    @test cfg.spatial_stride == 1
    @test cfg.derive_products == true
    @test cfg.memory_gb === nothing
    @test cfg.lambda_amp == 100.0
    @test cfg.lambda_mu == 100.0
    @test cfg.sig_init == 4.0
    @test isfile(cfg.config_path)
    @test wrapper.upstream_rohsa_source_dir() == normpath(joinpath(@__DIR__, "..", "rohsa", "upstream", "ROHSA"))
    @test wrapper.upstream_rohsa_binary() == joinpath(wrapper.upstream_rohsa_source_dir(), "src", "ROHSA")

    override_cfg, no_adjoint_cfg = mktempdir() do dir
        config_text = """
        [input]
        input = "data/exp_raw/component.fits"
        output_dir = "validation_output/rohsa_example"
        crop_size = 64
        crop_x0 = 1
        crop_y0 = 1
        reuse = true
        write_visuals = false

        [rohsa]
        n_gauss = [4, 6]
        lambda_amp = 250.0
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
        """
        load(text) = begin
            config_path = joinpath(dir, "rohsa.toml")
            write(config_path, text)
            wrapper.parse_args(["--config", config_path])
        end
        load(config_text), load(config_text * "regul_grad_adjoint = false\n")
    end
    @test override_cfg.n_gauss == [4, 6]
    @test override_cfg.lambda_amp == 250.0
    @test override_cfg.lambda_mu == 100.0
    @test override_cfg.sig_init == 4.0
    @test override_cfg.lb_sig == 1.0
    @test override_cfg.write_visuals == false
    @test override_cfg.spatial_stride == 1
    @test override_cfg.derive_products == true
    # One run per listed N, so both fit and neither overwrites the other.
    candidates = wrapper.resolved_rohsa_candidates(override_cfg)
    @test [c.n_gauss for c in candidates] == [4, 6]
    @test all(isempty(c.id) for c in candidates)
    @test_throws ArgumentError wrapper.parse_args(["--unknown-option", "value"])

    laptop_cfg = wrapper.parse_args(["--config", joinpath(@__DIR__, "..", "configs", "rohsa", "full_cube_laptop_stride8.toml")])
    @test laptop_cfg.spatial_stride == 8
    @test laptop_cfg.derive_products == false
    @test laptop_cfg.write_region_fits == false
    estimate = wrapper.rohsa_memory_estimate(453, 453, 106, 6, 10)
    @test estimate.padded_side == 512
    @test estimate.peak_gib < 10

    mktempdir() do dir
        input_path = joinpath(dir, "input.fits")
        cube = reshape(Float32.(1:48), 4, 4, 3)
        header = FITSHeader(
            [
                "SIMPLE", "BITPIX", "NAXIS", "NAXIS1", "NAXIS2", "NAXIS3",
                "CRPIX1", "CRPIX2", "CRPIX3",
                "CRVAL1", "CRVAL2", "CRVAL3",
                "CDELT1", "CDELT2", "CDELT3",
            ],
            [
                true, -32, 3, 4, 4, 3,
                1.0, 1.0, 1.0,
                0.0, 0.0, 0.0,
                1.0, 1.0, 1.0,
            ],
            fill("", 15),
        )
        f = FITS(input_path, "w")
        try
            write(f, cube, header=header)
        finally
            close(f)
        end
        region = wrapper.InputRegion(
            label="full_stride2",
            x0=1,
            y0=1,
            nx=2,
            ny=2,
            source_nx=4,
            source_ny=4,
            x_stride=2,
            y_stride=2,
        )
        products = wrapper.write_region_products(input_path, header, region, joinpath(dir, "out");
                                                 retain_cube=false, write_region_fits=false)
        @test products.cube === nothing
        @test !isfile(products.fits_path)
        @test isfile(products.rohsa_input_path)
        @test isfile(products.rms_path)
        @test first(readlines(products.rohsa_input_path)) == "3\t2\t2"
        @test_throws ArgumentError wrapper.write_region_products(input_path, header, region, joinpath(dir, "out_k");
                                                                 retain_cube=false, write_region_fits=false,
                                                                 spectral_kernel="hanning")
    end

    # spectral_kernel: the in-memory writer smooths, records the kernel in the crop
    # header, and scales the raw noise estimate by sqrt(sum k^2).
    mktempdir() do dir
        input_path = joinpath(dir, "noise.fits")
        noise = Float32.(0.3 .* randn(Random.MersenneTwister(3), 16, 16, 60))
        header = FITSHeader(
            ["SIMPLE", "BITPIX", "NAXIS", "NAXIS1", "NAXIS2", "NAXIS3",
             "CRPIX1", "CRPIX2", "CRPIX3", "CRVAL1", "CRVAL2", "CRVAL3", "CDELT1", "CDELT2", "CDELT3"],
            [true, -32, 3, 16, 16, 60, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            fill("", 15),
        )
        f = FITS(input_path, "w")
        try
            write(f, noise, header=header)
        finally
            close(f)
        end
        region = wrapper.InputRegion(label="full", x0=1, y0=1, nx=16, ny=16)
        plain = wrapper.write_region_products(input_path, header, region, joinpath(dir, "plain"))
        hann = wrapper.write_region_products(input_path, header, region, joinpath(dir, "hann");
                                             spectral_kernel="hanning")
        @test hann.sigma ≈ plain.sigma * sqrt(0.375)
        @test hann.kernel_sigma_chan ≈ sqrt(0.5)
        @test read_header(FITS(hann.fits_path)[1])["SPKERNEL"] == "hanning"

        # holdout_mask: held-out pixels get noise 0 (ROHSA skips their data term), the rest keep the cube's level.
        mask = zeros(Float32, 16, 16); mask[3:4, 5:8] .= 1  # mask[x, y]
        mask_path = joinpath(dir, "holdout.fits")
        FITS(mask_path, "w") do f
            write(f, mask)
        end
        held = wrapper.write_region_products(input_path, header, region, joinpath(dir, "held"); holdout_mask=mask_path)
        @test held.sigma ≈ plain.sigma
        vals = Dict((parse(Int, a), parse(Int, b)) => parse(Float64, c)
                    for (a, b, c) in split.(readlines(held.rms_path)[2:end], '\t'))
        @test count(==(0.0), values(vals)) == 8
        @test vals[(4, 2)] == 0.0          # (y, x) 0-based = mask[3, 5]
        @test vals[(0, 0)] ≈ plain.sigma
    end

    rows = [
        Dict{String, Any}("n_gauss" => 4, "total_model_over_data" => 1.2,
                          "moment_relative_rmse" => 0.20, "moment_mean_relative_bias" => 0.10,
                          "moment_residual_skew" => 1.0, "selection_score" => 0.60,
                          "product_root" => "n4", "run_trace" => "run4.toml",
                          "products_trace" => "products4.toml"),
        Dict{String, Any}("n_gauss" => 6, "total_model_over_data" => 1.01,
                          "moment_relative_rmse" => 0.03, "moment_mean_relative_bias" => 0.01,
                          "moment_residual_skew" => 0.1, "selection_score" => 0.06,
                          "product_root" => "n6", "run_trace" => "run6.toml",
                          "products_trace" => "products6.toml"),
    ]
    @test wrapper.select_best_rohsa_n(rows)["n_gauss"] == 6
    selected_rows = wrapper.selected_rohsa_summary_rows(Dict("full" => wrapper.select_best_rohsa_n(rows)))
    @test only(selected_rows)["region"] == "full"
    @test only(selected_rows)["selected_n_gauss"] == 6
    @test only(selected_rows)["product_root"] == "n6"

    mktempdir() do dir
        params_path = joinpath(dir, "parameters.txt")
        dat_path = joinpath(dir, "cube.dat")
        rms_path = joinpath(dir, "rms.dat")
        out_path = joinpath(dir, "out.dat")
        timeout_path = joinpath(dir, "timeout.dat")
        cand = only(filter(c -> c.n_gauss == 6, candidates))
        wrapper.write_rohsa_parameters(params_path, dat_path, rms_path, out_path, timeout_path, override_cfg, cand)
        text = read(params_path, String)
        @test occursin("lambda_amp = 250.0d0", text)
        @test occursin("lambda_mu = 100.0d0", text)
        @test occursin("sig_init = 4.0d0", text)
        @test occursin("noise = .false.", text)
        # On by default. "false" writes nothing, so ROHSA binaries without the switch still read the file.
        @test occursin("regul_grad_adjoint = .true.", text)
        wrapper.write_rohsa_parameters(params_path, dat_path, rms_path, out_path, timeout_path, no_adjoint_cfg, cand)
        @test !occursin("regul_grad_adjoint", read(params_path, String))
    end

    # Per-pixel noise map for ROHSA: header "ny nx", rows "y x sigma" (0-based), bad pixels written as 0
    # (ROHSA then skips their data term), Julia sigma[x, y] order.
    mktempdir() do dir
        sigma = [0.2 NaN; -1.0 0.3]  # sigma[x, y]
        path = wrapper.write_rms_map(joinpath(dir, "rms_map.dat"), sigma)
        lines = readlines(path)
        @test lines[1] == "2\t2"
        vals = Dict((parse(Int, a), parse(Int, b)) => parse(Float64, c) for (a, b, c) in split.(lines[2:end], '\t'))
        @test vals[(0, 0)] == 0.2  # (y, x) = (0, 0) -> sigma[1, 1]
        @test vals[(1, 0)] == 0.0  # sigma[1, 2] = NaN
        @test vals[(0, 1)] == 0.0  # sigma[2, 1] = -1
        @test vals[(1, 1)] == 0.3
    end

    # Binary ROHSA input contract with mod_inout.f90 read_cube: int32 nv, ny,
    # nx, then Float64 in Fortran cube(v, y, x) order.
    mktempdir() do dir
        cube = reshape(Float32.(1:24), 2, 3, 4)  # Julia (x, y, v)
        path = wrapper.write_cube_bin(joinpath(dir, "crop.bin"), cube)
        open(path) do io
            @test [read(io, Int32) for _ in 1:3] == Int32[4, 3, 2]
            vals = reinterpret(Float64, read(io))
            @test vals == vec(Float64.(permutedims(cube, (3, 2, 1))))
        end
    end
end
