using Astronomy
using FITSIO
using Test

@testset "ROHSA backend derives physical component products" begin
    header = FITSHeader(
        [
            "SIMPLE", "BITPIX", "NAXIS", "NAXIS1", "NAXIS2", "NAXIS3",
            "CRPIX1", "CRPIX2", "CRPIX3",
            "CRVAL1", "CRVAL2", "CRVAL3",
            "CDELT1", "CDELT2", "CDELT3",
            "CUNIT3",
        ],
        [
            true, -32, 3, 2, 2, 3,
            1.0, 1.0, 1.0,
            0.0, 0.0, 0.0,
            1.0, 1.0, 1000.0,
            "m/s",
        ],
        fill("", 16),
    )
    velocities, dv = Astronomy.velocity_axis_kms(header, 3)
    @test velocities == [-1.0, 0.0, 1.0]
    @test dv == 1.0
    @test Astronomy.doppler_temperature_k(2.0) ≈ 484.4

    params = zeros(Float64, 6, 2, 2)
    params[1, :, :] .= 4.0
    params[2, :, :] .= 1.0
    params[3, :, :] .= 1.0
    params[4, :, :] .= 2.0
    params[5, :, :] .= 1.0
    params[6, :, :] .= 8.0
    components = Astronomy.reconstruct_rohsa_components(params, 3)
    data_cube = components[1] .+ components[2]

    rows = Astronomy.rohsa_component_rows("region", 2, params, components, velocities, dv)
    @test rows[1]["width_group"] == "narrow"
    @test rows[1]["td_group"] == "td_lt_500k"
    @test rows[2]["width_group"] == "broad"
    @test rows[2]["td_group"] == "td_ge_5000k"
    @test rows[1]["weighted_sigma_kms"] ≈ 1.0
    @test rows[2]["weighted_t_doppler_k"] ≈ 121.1 * 64.0

    mktempdir() do dir
        products = Astronomy.write_rohsa_decomposition_products(
            dir,
            "region",
            header,
            params,
            components,
            data_cube;
            n_gauss=2,
        )
        @test isfile(products.model_path)
        @test isfile(products.cube_paths["model"])
        @test isfile(products.cube_paths["residual"])
        @test isfile(products.cube_paths["component_g1"])
        @test isfile(products.group_paths["narrow"])
        @test isfile(products.group_paths["td_lt_500k"])
        @test isfile(products.map_paths["component_g1_t_doppler_k"])
        @test isfile(products.map_paths["td_lt_500k_fraction"])
        @test products.fit_row["total_model_over_data"] ≈ 1.0
        @test products.group_nhi["td_lt_500k_fraction"] == ones(2, 2) .* (products.group_nhi["td_lt_500k"][1, 1] / products.model_nhi[1, 1])

        labeled = Astronomy.write_rohsa_decomposition_products(
            dir,
            "region",
            header,
            params,
            components,
            data_cube;
            n_gauss=2,
            run_label="rohsa_c001_n2",
        )
        @test occursin("rohsa_c001_n2", labeled.root)
        @test occursin("rohsa_c001_n2_model.fits", labeled.model_path)

        if Sys.which("fpack") === nothing
            @test_throws ArgumentError Astronomy.write_rohsa_decomposition_products(
                dir,
                "region",
                header,
                params,
                components,
                data_cube;
                n_gauss=2,
                run_label="rohsa_c002_n2",
                fits_extension=".fits.fz",
            )
        else
            compressed = Astronomy.write_rohsa_decomposition_products(
                dir,
                "region",
                header,
                params,
                components,
                data_cube;
                n_gauss=2,
                run_label="rohsa_c002_n2",
                fits_extension=".fits.fz",
            )
            @test endswith(compressed.model_path, ".fits.fz")
            @test isfile(compressed.model_path)
            @test Astronomy.read_fits_image(compressed.model_path) == compressed.model
        end
    end
end

@testset "ROHSA output reader skips comments and validates row count" begin
    mktempdir() do dir
        path = joinpath(dir, "rohsa.dat")
        open(path, "w") do io
            println(io, "# y x amp mu sig")
            println(io, "0 0 1.5 2.0 3.0")
        end
        params = Astronomy.read_rohsa_parameters(path, 1, 1, 1)
        @test params[:, 1, 1] == [1.5, 1.0, 3.0]
        @test_throws ArgumentError Astronomy.read_rohsa_parameters(path, 2, 1, 1)
    end
end
