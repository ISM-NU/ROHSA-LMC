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

    override_cfg, no_adjoint_cfg, budget_cfg, noise_cfg = mktempdir() do dir
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
        @test_throws ArgumentError load(config_text * "checkpoint_at = [0]\n")
        load(config_text), load(config_text * "regul_grad_adjoint = false\n"),
        load(config_text * "maxiter_last = 5000\ncheckpoint_at = [800, 100]\n"),
        load(config_text * "noise_estimate = \"line_free\"\nnoise_per_pixel = true\nnoise_channels = \"1-10,95-106\"\n" *
             "noise_floor_quantile = 0.01\nmissing_data = \"skip\"\n")
    end
    @test override_cfg.n_gauss == [4, 6]
    @test (override_cfg.noise_estimate, override_cfg.noise_per_pixel, override_cfg.noise_channels) ==
          ("channel_difference", false, "")
    @test (noise_cfg.noise_estimate, noise_cfg.noise_per_pixel, noise_cfg.noise_channels) ==
          ("line_free", true, "1-10,95-106")
    @test (override_cfg.missing_data, override_cfg.noise_floor_quantile) == ("zero", 0.0)
    @test (noise_cfg.missing_data, noise_cfg.noise_floor_quantile) == ("skip", 0.01)
    @test override_cfg.lambda_amp == 250.0
    @test override_cfg.lambda_mu == 100.0
    @test override_cfg.sig_init == 4.0
    @test override_cfg.lb_sig == 1.0
    @test override_cfg.write_visuals == false
    @test override_cfg.spatial_stride == 1
    @test override_cfg.derive_products == true
    # A plain n_gauss list is an error; several N need the [auto_config] n_gauss axis.
    @test_throws ArgumentError wrapper.resolved_rohsa_candidates(override_cfg)
    single = only(wrapper.resolved_rohsa_candidates(cfg))
    @test single.n_gauss == 6 && isempty(single.id)
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

        # Default output unchanged by the new options: the scalar from the channel-difference estimate, written by the
        # scalar writer byte for byte.
        @test plain.sigma == wrapper.robust_channel_difference_noise(noise)
        ref_rms = joinpath(dir, "ref_rms.dat")
        wrapper.write_rms_map(ref_rms, 16, 16, plain.sigma)
        @test read(ref_rms) == read(plain.rms_path)

        # noise_estimate = line_free: scalar by default, the per-pixel map with noise_per_pixel; held-out pixels
        # still get 0; bad option values are rejected.
        lf = wrapper.write_region_products(input_path, header, region, joinpath(dir, "lf"); noise_estimate="line_free")
        @test isapprox(lf.sigma, 0.3; rtol=0.1)
        lfvals = unique(parse(Float64, split(l, '\t')[3]) for l in readlines(lf.rms_path)[2:end])
        @test length(lfvals) == 1
        lfp = wrapper.write_region_products(input_path, header, region, joinpath(dir, "lfp");
                                            noise_estimate="line_free", noise_per_pixel=true, holdout_mask=mask_path)
        pvals = [parse(Float64, split(l, '\t')[3]) for l in readlines(lfp.rms_path)[2:end]]
        @test count(==(0.0), pvals) == 8
        @test length(unique(filter(>(0), pvals))) > 1
        @test all(v -> v == 0 || isapprox(v, 0.3; rtol=0.25), pvals)
        @test_throws ArgumentError wrapper.write_region_products(input_path, header, region, joinpath(dir, "bad");
                                                                 noise_estimate="lstd")
        @test_throws ArgumentError wrapper.write_region_products(input_path, header, region, joinpath(dir, "bad2");
                                                                 noise_per_pixel=true)
        @test plain.range_path == "" && lfp.range_path == ""   # no range file unless missing_data = skip

        # missing_data = skip: live channel ranges from the raw cube go to live_channel_range.dat, pixels without data
        # get noise 0 in rms_map.dat (scalar or per-pixel noise alike); the line-free map can be floored.
        holes = copy(noise)
        holes[1:3, 1:3, :] .= NaN32          # beyond the survey edge
        holes[9:end, :, 1:20] .= 0.0f0       # an observing block without channels 1-20
        holes_path = joinpath(dir, "holes.fits")
        FITS(holes_path, "w") do f
            write(f, holes, header=header)
        end
        sk = wrapper.write_region_products(holes_path, header, region, joinpath(dir, "skip"); missing_data="skip")
        @test isfile(sk.range_path) && basename(sk.range_path) == "live_channel_range.dat"
        rows = split.(readlines(sk.range_path)[2:end], '\t')
        rng_of = Dict((parse(Int, r[1]), parse(Int, r[2])) => (parse(Int, r[3]), parse(Int, r[4])) for r in rows)
        @test rng_of[(0, 0)] == (0, -1)            # (y, x) 0-based: holes[1, 1] has no data
        @test rng_of[(0, 8)] == (20, 59)           # holes[9, 1]: channels 21-60, 0-based 20-59
        @test rng_of[(5, 5)] == (0, 59)
        skv = Dict((parse(Int, a), parse(Int, b)) => parse(Float64, c)
                   for (a, b, c) in split.(readlines(sk.rms_path)[2:end], '\t'))
        @test count(==(0.0), values(skv)) == 9 && skv[(2, 2)] == 0.0
        @test all(v -> v == 0 || v == sk.sigma, values(skv))
        skl = wrapper.write_region_products(holes_path, header, region, joinpath(dir, "skipl"); missing_data="skip",
                                            noise_estimate="line_free", noise_per_pixel=true, noise_floor_quantile=0.25)
        lv = [parse(Float64, split(l, '\t')[3]) for l in readlines(skl.rms_path)[2:end]]
        @test count(==(0.0), lv) == 9
        live = filter(>(0), lv)
        @test count(==(minimum(live)), live) >= length(live) ÷ 5   # about a quarter of the pixels sit on the floor
        @test_throws ArgumentError wrapper.write_region_products(input_path, header, region, joinpath(dir, "bad3");
                                                                 missing_data="nan")
        @test_throws ArgumentError wrapper.write_region_products(input_path, header, region, joinpath(dir, "bad4");
                                                                 noise_floor_quantile=1.0)
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
        cand = wrapper.RohsaCandidate("", 6, override_cfg.lambda_amp, override_cfg.lambda_mu, override_cfg.lambda_sig,
                                      override_cfg.lambda_var_amp, override_cfg.lambda_var_mu, override_cfg.lambda_var_sig)
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
        # Budget-study switches: absent unless set; checkpoints go next to the output, sorted.
        @test !occursin("maxiter_last", text) && !occursin("checkpoint", text)
        @test budget_cfg.maxiter_last == 5000 && budget_cfg.checkpoint_at == [100, 800]
        wrapper.write_rohsa_parameters(params_path, dat_path, rms_path, out_path, timeout_path, budget_cfg, cand)
        text = read(params_path, String)
        @test occursin("maxiter_last = 5000", text)
        @test occursin("checkpoint_at = 100,800", text)
        @test occursin("checkpoint_prefix = '$(joinpath(dir, "checkpoints", "rohsa_n6"))'", text)
        # Live channel range file: emitted only when given.
        @test !occursin("filename_range", text)
        wrapper.write_rohsa_parameters(params_path, dat_path, rms_path, out_path, timeout_path, override_cfg, cand;
                                       range_path=joinpath(dir, "live_channel_range.dat"))
        @test occursin("filename_range = '$(joinpath(dir, "live_channel_range.dat"))'", read(params_path, String))
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
