using Astronomy
using RohsaLMC
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
    # CRPIX3 = 1, CRVAL3 = 0: the first channel is at 0 km/s (FITS convention).
    @test velocities == [0.0, 1.0, 2.0]
    @test dv == 1.0
    @test RohsaLMC.doppler_temperature_k(2.0) ≈ 484.4

    params = zeros(Float64, 6, 2, 2)
    params[1, :, :] .= 4.0
    params[2, :, :] .= 1.0
    params[3, :, :] .= 1.0
    params[4, :, :] .= 2.0
    params[5, :, :] .= 1.0
    params[6, :, :] .= 8.0
    components = RohsaLMC.reconstruct_rohsa_components(params, 3)
    data_cube = components[1] .+ components[2]

    rows = RohsaLMC.rohsa_component_rows("region", 2, params, components, velocities, dv)
    @test rows[1]["width_group"] == "narrow"
    @test rows[1]["td_group"] == "td_lt_500k"
    @test rows[2]["width_group"] == "broad"
    @test rows[2]["td_group"] == "td_ge_5000k"
    @test rows[1]["weighted_sigma_kms"] ≈ 1.0
    # params hold zero-based centres: mu = 1 is the second channel (1 km/s here).
    @test rows[1]["weighted_mu_kms"] ≈ 1.0
    @test rows[2]["weighted_t_doppler_k"] ≈ 121.1 * 64.0

    mktempdir() do dir
        products = RohsaLMC.write_rohsa_decomposition_products(
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

        # Grouping is per pixel: a component narrow in one pixel and broad in the
        # others counts as cold only where it is narrow.
        mixed = copy(params)
        mixed[6, :, :] .= 8.0
        mixed[6, 1, 1] = 1.0
        mixed_components = RohsaLMC.reconstruct_rohsa_components(mixed, 3)
        mixed_products = RohsaLMC.write_rohsa_decomposition_products(
            dir, "region", header, mixed, mixed_components, mixed_components[1] .+ mixed_components[2];
            n_gauss=2, write_cubes=false, run_label="rohsa_mixed",
        )
        @test mixed_products.group_nhi["td_lt_500k_fraction"][1, 1] ≈ 1.0
        @test mixed_products.group_nhi["td_lt_500k_fraction"][2, 2] ≈
              mixed_products.group_nhi["td_lt_500k"][2, 2] / mixed_products.model_nhi[2, 2]
        @test mixed_products.group_nhi["td_lt_500k_fraction"][2, 2] < 1.0

        # Cubes off (the CLI default): maps and metrics only, no 3D files.
        mktempdir() do slim_dir
            slim = RohsaLMC.write_rohsa_decomposition_products(
                slim_dir, "region", header, params, components, data_cube;
                n_gauss=2, write_cubes=false,
            )
            @test isempty(slim.cube_paths)
            @test slim.model_path == ""
            @test all(isempty, values(slim.component_paths))
            @test isfile(slim.map_paths["component_g1_nhi_cm2"])
            @test isfile(slim.map_paths["model_nhi_cm2"])
            @test isfile(slim.map_paths["td_lt_500k_fraction"])
            @test slim.fit_row["total_model_over_data"] ≈ 1.0
            @test isempty(filter(f -> endswith(f, ".fits") && occursin("_group", f), readdir(slim_dir)))
        end

        labeled = RohsaLMC.write_rohsa_decomposition_products(
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
            @test_throws ArgumentError RohsaLMC.write_rohsa_decomposition_products(
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
            compressed = RohsaLMC.write_rohsa_decomposition_products(
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
        params = RohsaLMC.read_rohsa_parameters(path, 1, 1, 1)
        @test params[:, 1, 1] == [1.5, 1.0, 3.0]
        # End to end: ROHSA's mu is 1-based, so mu = 2.0 in the file is the second
        # channel, at 1 km/s on an axis with CRVAL3 = 0, CRPIX3 = 1, CDELT3 = 1 km/s.
        header = FITSHeader(["CRPIX3", "CRVAL3", "CDELT3", "CUNIT3"], [1.0, 0.0, 1000.0, "m/s"], fill("", 4))
        velocities, dv = Astronomy.velocity_axis_kms(header, 4)
        components = RohsaLMC.reconstruct_rohsa_components(params, 4)
        rows = RohsaLMC.rohsa_component_rows("r", 1, params, components, velocities, dv)
        @test rows[1]["weighted_mu_kms"] ≈ 1.0
        @test argmax(vec(components[1][1, 1, :])) == 2
        @test_throws ArgumentError RohsaLMC.read_rohsa_parameters(path, 2, 1, 1)
    end
end

@testset "ROHSA map products state their own BUNIT" begin
    # Cube header with BUNIT = K, as the H I cubes have: maps must not inherit it.
    header = FITSHeader(
        ["SIMPLE", "BITPIX", "NAXIS", "NAXIS1", "NAXIS2", "NAXIS3", "CRPIX3", "CRVAL3", "CDELT3", "CUNIT3", "BUNIT"],
        [true, -32, 3, 2, 2, 3, 1.0, 0.0, 1000.0, "m/s", "K"],
        fill("", 11),
    )
    params = zeros(Float64, 3, 2, 2)
    params[1, :, :] .= 4.0
    params[2, :, :] .= 1.0
    params[3, :, :] .= 1.0
    components = RohsaLMC.reconstruct_rohsa_components(params, 3)
    data_cube = components[1]
    bunit(path) = FITS(path) do f
        h = read_header(f[1])
        haskey(h, "BUNIT") ? strip(h["BUNIT"]) : nothing  # CFITSIO stores "" as " "
    end
    expected = Dict(
        "component_g1_amp_k" => "K",
        "component_g1_mu_kms" => "km/s",
        "component_g1_sigma_kms" => "km/s",
        "component_g1_t_doppler_k" => "K",
        "component_g1_nhi_cm2" => "cm-2",
        "narrow_nhi_cm2" => "cm-2",
        "td_lt_500k_nhi_cm2" => "cm-2",
        "model_nhi_cm2" => "cm-2",
        "data_nhi_cm2" => "cm-2",
        "residual_nhi_cm2" => "cm-2",
        "residual_moment0_kkms" => "K km/s",
        "narrow_fraction" => "",
    )
    mktempdir() do dir
        products = RohsaLMC.write_rohsa_decomposition_products(dir, "region", header, params, components, data_cube; n_gauss=1)
        for (name, unit) in expected
            @test bunit(products.map_paths[name]) == unit
        end
        # Every written map is covered, so a new product cannot slip in with the cube's K.
        @test Set(keys(products.map_paths)) ⊇ Set(keys(expected))
        for (name, path) in products.map_paths
            @test RohsaLMC.rohsa_map_bunit(name) !== nothing
            @test bunit(path) == RohsaLMC.rohsa_map_bunit(name)
        end
        # The cubes are brightness and keep the cube's unit.
        @test bunit(products.model_path) == "K"
        @test bunit(products.cube_paths["residual"]) == "K"
    end

    # A product with no known unit gets no BUNIT, not the cube's K.
    @test RohsaLMC.rohsa_map_bunit("something_new") === nothing
    h = RohsaLMC.header_with_bunit(header, nothing)
    @test !haskey(h, "BUNIT")
    @test header["BUNIT"] == "K"
end
