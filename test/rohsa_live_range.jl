using Test
using Random
using RohsaLMC: read_rohsa_parameters

# Contract of the live channel range (mod_array.f90 filename_range, used by missing_data = "skip"):
# 1. Without a range file nothing changes, and a range file that marks every channel live gives the same output byte
#    for byte (the change is a no-op on full ranges).
# 2. Channels outside a pixel's range carry no data term: on a cube where a block of channels is zeroed in half of the
#    field (as in a mosaic whose pointings cover different channel ranges), the fit with ranges recovers the injected
#    Gaussians about as well as the fit of the cube without the zeroed block, and better than a fit that treats the
#    zeros as data.
# 3. A pixel without data is fitted from the smoothing term alone: its amplitudes are not pinned to 0.
@testset "ROHSA live channel range" begin
    if Sys.which("gfortran") === nothing || Sys.which("make") === nothing
        @test_skip "gfortran or make not available"
    else
        src = joinpath(@__DIR__, "..", "rohsa", "upstream", "ROHSA", "src")
        mktempdir() do dir
            build = joinpath(dir, "src")
            mkpath(build)
            for f in readdir(src)
                (endswith(f, ".f") || endswith(f, ".f90") || endswith(f, ".h") || f == "Makefile") &&
                    cp(joinpath(src, f), joinpath(build, f))
            end
            run(pipeline(Cmd(`make all`; dir=build); stdout=devnull, stderr=devnull))
            bin = joinpath(build, "ROHSA")

            nv, ny, nx, ng = 80, 32, 32, 3
            gauss(v, a, mu, s) = a * exp(-(v - mu)^2 / (2s^2))
            truth = Array{Float64}(undef, nv, ny, nx)   # Fortran order (v, y, x)
            for x in 1:nx, y in 1:ny, v in 1:nv
                truth[v, y, x] = gauss(v, 3.0, 30.0 + 0.15x, 3.0) + gauss(v, 1.5, 50.0 - 0.1y, 2.0)
            end
            rng = MersenneTwister(12)
            sigma = 0.1
            observed = truth .+ sigma .* randn(rng, nv, ny, nx)

            function write_cube(name, c)
                path = joinpath(dir, name)
                open(path, "w") do io
                    write(io, Int32(nv), Int32(ny), Int32(nx))
                    write(io, c)
                end
                return path
            end
            function write_map(name, m)          # m[y, x]
                path = joinpath(dir, name)
                open(path, "w") do io
                    println(io, "$ny\t$nx")
                    for y in 1:ny, x in 1:nx
                        println(io, "$(y - 1)\t$(x - 1)\t$(m[y, x])")
                    end
                end
                return path
            end
            function write_ranges(name, lo, hi)  # 1-based; written 0-based like the wrapper
                path = joinpath(dir, name)
                open(path, "w") do io
                    println(io, "$ny\t$nx")
                    for y in 1:ny, x in 1:nx
                        println(io, "$(y - 1)\t$(x - 1)\t$(lo[y, x] - 1)\t$(hi[y, x] - 1)")
                    end
                end
                return path
            end
            function fit(name, cube_path, rms_path; noise=".true.", range_path="", maxiter=200)
                out = joinpath(dir, name)
                mkpath(out)
                params = joinpath(out, "parameters.txt")
                extra = isempty(range_path) ? "" : ",filename_range = '$range_path'"
                write(params, """
                &user_parameters
                    filename = '$cube_path'
                    ,fileout = '$(joinpath(out, "out.dat"))'
                    ,timeout = '$(joinpath(out, "timestep.dat"))'
                    ,filename_noise = '$rms_path'
                    ,n_gauss = $ng
                    ,lambda_amp = 10.0d0
                    ,lambda_mu = 10.0d0
                    ,lambda_sig = 10.0d0
                    ,lambda_var_sig = 0.0d0
                    ,amp_fact_init = 0.66d0
                    ,sig_init = 4.0d0
                    ,lb_sig_init = 1.0d0
                    ,ub_sig_init = 12.0d0
                    ,lb_sig = 1.0d0
                    ,ub_sig = 100.0d0
                    ,init_option = 'mean'
                    ,maxiter_init = 15000
                    ,maxiter = $maxiter
                    ,m = 10
                    ,noise = $noise
                    ,regul = .true.
                    ,descent = .true.
                    ,lstd = 1
                    ,ustd = 10
                    ,iprint = -1
                    ,iprint_init = -1
                    ,save_grid = .false.
                    $extra
                    /
                """)
                cmd = addenv(Cmd(`$bin $params`; dir=out), "OMP_NUM_THREADS" => "2")
                run(pipeline(cmd; stdout=joinpath(out, "stdout.log"), stderr=joinpath(out, "stderr.log")))
                return joinpath(out, "out.dat")
            end
            # Model cube (v, y, x) from an output file; read_rohsa_parameters gives (3N, x, y) with 0-based centres.
            function model(out)
                p = read_rohsa_parameters(out, ng, nx, ny)
                m = zeros(nv, ny, nx)
                for x in 1:nx, y in 1:ny, g in 1:ng, v in 1:nv
                    m[v, y, x] += gauss(v - 1, p[3g - 2, x, y], p[3g - 1, x, y], max(p[3g, x, y], 1e-6))
                end
                return m, p
            end

            # 1. No range file vs a range file with every channel live: identical output (scalar and per-pixel noise).
            full_cube = write_cube("full.bin", observed)
            flat = write_map("flat.dat", fill(sigma, ny, nx))
            varied = write_map("varied.dat", [sigma * (1 + 0.3 * sin(x / 5) * cos(y / 7)) for y in 1:ny, x in 1:nx])
            all_live = write_ranges("all_live.dat", ones(Int, ny, nx), fill(nv, ny, nx))
            for (tag, rms, noise) in (("flat", flat, ".true."), ("varied", varied, ".true."), ("own", flat, ".false."))
                a = fit("plain_$tag", full_cube, rms; noise=noise, maxiter=40)
                b = fit("alllive_$tag", full_cube, rms; noise=noise, range_path=all_live, maxiter=40)
                @test read(a) == read(b)
            end

            # 2. Zero channels 1-28 in the right half (they cut the lower wing of the first line), and remove a 6x6
            #    corner entirely (no data: zeros, noise 0, empty range).
            zeroed = copy(observed)
            zeroed[1:28, :, 17:end] .= 0
            nodata = falses(ny, nx); nodata[1:6, 1:6] .= true
            for x in 1:nx, y in 1:ny
                nodata[y, x] && (zeroed[:, y, x] .= 0)
            end
            lo = [nodata[y, x] ? 1 : (x >= 17 ? 29 : 1) for y in 1:ny, x in 1:nx]
            hi = [nodata[y, x] ? 0 : nv for y in 1:ny, x in 1:nx]
            zcube = write_cube("zeroed.bin", zeroed)
            rms0 = write_map("rms_nodata.dat", [nodata[y, x] ? 0.0 : sigma for y in 1:ny, x in 1:nx])
            rng_path = write_ranges("ranges.dat", lo, hi)

            m_full, _ = model(fit("full", full_cube, flat))
            m_zero, p_zero = model(fit("zeros_as_data", zcube, rms0))
            m_rng, p_rng = model(fit("ranges", zcube, rms0; range_path=rng_path))

            # Error against the injected (noise-free) model over all channels, in the zeroed half.
            cut = [!nodata[y, x] && x >= 17 for y in 1:ny, x in 1:nx]
            err(m) = sqrt(sum(abs2, (m .- truth)[:, cut]) / (nv * count(cut)))
            e_full, e_zero, e_rng = err(m_full), err(m_zero), err(m_rng)
            println("live range test: truth RMSE (K, noise 0.1) in the zeroed half: full cube $e_full, ",
                    "zeros as data $e_zero, with ranges $e_rng")
            # 2026-11-22: 0.0275 / 0.0600 / 0.0296
            @test e_rng <= 1.25 * e_full
            @test e_zero > 1.5 * e_rng
            # Column of the model in the zeroed half (the zeroed block cuts the first line's lower wing).
            col(m) = sum(m[:, cut]) / sum(truth[:, cut]) - 1
            println("live range test: column error: full $(col(m_full)), zeros as data $(col(m_zero)), ranges $(col(m_rng))")
            # 2026-11-22: -0.002 / -0.032 / -0.004
            @test abs(col(m_rng)) <= abs(col(m_full)) + 0.005
            @test abs(col(m_zero)) > 2 * abs(col(m_rng))

            # 3. Pixels without data: pinned at 0 when their zeros are data; with ranges the smoothing term fills them
            #    with a model close to their neighbours'.
            nd = [nodata[y, x] for y in 1:ny, x in 1:nx]
            @test maximum(abs, m_zero[:, nd]) < 0.05
            ratio = sum(m_rng[:, nd]) / sum(truth[:, nd])
            println("live range test: model column / truth at pixels without data, with ranges: $ratio")  # 1.006
            @test 0.5 < ratio < 1.5
        end
    end
end
