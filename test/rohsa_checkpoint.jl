using Test
using Random
using RohsaLMC: read_rohsa_parameters, read_rohsa_checkpoint

# Contract of the budget-study switches (mod_array.f90 maxiter_last, checkpoint_at): the evaluation limit is only
# checked by the driver and never enters L-BFGS-B, so the iterates do not depend on it. A checkpoint written at K
# evaluations of the full-resolution level is therefore exactly the fit a run with limit K returns. Budget curves
# are read off one long run on this guarantee.
@testset "ROHSA checkpoint at K equals a run stopped at K" begin
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

            nv, ny, nx = 40, 16, 16
            rng = MersenneTwister(2)
            cube = Array{Float64}(undef, nv, ny, nx)
            for x in 1:nx, y in 1:ny, v in 1:nv
                s = 0.0
                for (a, mu, sig) in ((4.0, 12.0 + 0.2x, 2.5), (2.0, 24.0 - 0.2y, 1.5))
                    s += a * exp(-(v - mu)^2 / (2sig^2))
                end
                cube[v, y, x] = s + 0.1randn(rng)
            end
            open(joinpath(dir, "cube.bin"), "w") do io
                write(io, Int32(nv), Int32(ny), Int32(nx))
                write(io, cube)
            end

            function fit(name, extra)
                out = joinpath(dir, name)
                mkpath(out)
                params = joinpath(out, "parameters.txt")
                write(params, """
                &user_parameters
                    filename = '$(joinpath(dir, "cube.bin"))'
                    ,fileout = '$(joinpath(out, "out.dat"))'
                    ,timeout = '$(joinpath(out, "timestep.dat"))'
                    ,n_gauss = 4
                    ,lambda_amp = 10.0d0
                    ,lambda_mu = 10.0d0
                    ,lambda_sig = 10.0d0
                    ,amp_fact_init = 0.66d0
                    ,sig_init = 4.0d0
                    ,lb_sig_init = 1.0d0
                    ,ub_sig_init = 12.0d0
                    ,lb_sig = 1.0d0
                    ,ub_sig = 100.0d0
                    ,init_option = 'mean'
                    ,maxiter_init = 15000
                    ,maxiter = 25
                    ,m = 10
                    ,noise = .false.
                    ,regul = .true.
                    ,descent = .true.
                    ,lstd = 1
                    ,ustd = 6
                    ,iprint = -1
                    ,iprint_init = -1
                    ,save_grid = .false.
                    $extra
                    /
                """)
                cmd = addenv(Cmd(`$(joinpath(build, "ROHSA")) $params`; dir=out), "OMP_NUM_THREADS" => "2")
                run(pipeline(cmd; stdout=joinpath(out, "stdout.log"), stderr=joinpath(out, "stderr.log")))
                return out
            end

            plain = fit("plain", "")
            long = fit("long", ",maxiter_last = 200\n,checkpoint_at = 10, 25\n,checkpoint_prefix = '$(joinpath(dir, "ck"))'")
            p25 = read_rohsa_parameters(joinpath(plain, "out.dat"), 4, nx, ny)
            c25, h = read_rohsa_checkpoint(joinpath(dir, "ck_e25.ckpt"))
            @test size(c25) == size(p25)
            @test h.evals >= 25 && h.n_gauss == 4 && length(h.b_params) == 4
            # out.dat is list-directed text, so allow its last-digit rounding.
            @test isapprox(c25, p25; rtol=1e-12, atol=1e-12)
            @test isfile(joinpath(dir, "ck_e10.ckpt"))
            # The long run went on: the final fit differs from the K = 25 one, and its log has trajectory lines.
            @test read_rohsa_parameters(joinpath(long, "out.dat"), 4, nx, ny) != p25
            log = read(joinpath(long, "stdout.log"), String)
            @test count("trajectory iter=", log) > 25
            @test occursin("checkpoint at=25", log)
            @test !occursin("trajectory", read(joinpath(plain, "stdout.log"), String))
        end
    end
end
