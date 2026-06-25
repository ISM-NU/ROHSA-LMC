function derive_rohsa_products_stage(region_dir, region_label::AbstractString, header, data_cube,
                                     rohsa_out::AbstractString, n_gauss::Int)
    params = Astronomy.read_rohsa_parameters(rohsa_out, n_gauss, size(data_cube, 1), size(data_cube, 2))
    components = Astronomy.reconstruct_rohsa_components(params, size(data_cube, 3))
    products = Astronomy.write_rohsa_decomposition_products(
        region_dir,
        region_label,
        header,
        params,
        components,
        data_cube;
        n_gauss=n_gauss,
    )
    trace_path = joinpath(region_dir, "rohsa_n$(n_gauss)_products_trace.toml")
    open(trace_path, "w") do io
        TOML.print(io, Dict(
            "stage" => "rohsa_product_derivation",
            "region" => region_label,
            "n_gauss" => n_gauss,
            "rohsa_output" => rohsa_out,
            "product_root" => products.root,
            "model_path" => products.model_path,
            "map_paths" => products.map_paths,
            "cube_paths" => products.cube_paths,
            "fit_metrics" => products.fit_row,
        ); sorted=true)
    end
    return (; products..., trace_path)
end
