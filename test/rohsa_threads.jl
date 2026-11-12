using Test
using Random

# Contract of the vendored Fortran: every OpenMP sum has a fixed order (omp_block.h,
# f_g_cube_fast), so the output does not depend on OMP_NUM_THREADS. Replicate runs
# therefore have to come from perturbed inputs, and a thread-count change must not
# change a result. n = 3 * 14 * 40^2 + 14 = 67214 is above the L-BFGS-B block length
# (65536), so the blocked sums are exercised too.
@testset "ROHSA output does not depend on OMP_NUM_THREADS" begin
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

            # Three spatially smooth Gaussians plus noise, cube(v, y, x) in Fortran order.
            nv, ny, nx = 48, 40, 40
            rng = MersenneTwister(1)
            cube = Array{Float64}(undef, nv, ny, nx)
            for x in 1:nx, y in 1:ny, v in 1:nv
                s = 0.0
                for (a, mu, sig) in ((5.0, 14.0 + 0.1x, 2.0), (3.0, 26.0 - 0.1y, 4.0), (1.5, 33.0, 1.5 + 0.02(x + y)))
                    s += a * exp(-(v - mu)^2 / (2sig^2))
                end
                cube[v, y, x] = s + 0.1randn(rng)
            end
            open(joinpath(dir, "cube.bin"), "w") do io
                write(io, Int32(nv), Int32(ny), Int32(nx))
                write(io, cube)
            end

            outputs = String[]
            for threads in (1, 3)
                out = joinpath(dir, "t$threads")
                mkpath(out)
                params = joinpath(out, "parameters.txt")
                write(params, """
                &user_parameters
                    filename = '$(joinpath(dir, "cube.bin"))'
                    ,fileout = '$(joinpath(out, "out.dat"))'
                    ,timeout = '$(joinpath(out, "timestep.dat"))'
                    ,n_gauss = 14
                    ,lambda_amp = 10.0d0
                    ,lambda_mu = 10.0d0
                    ,lambda_sig = 10.0d0
                    ,lambda_var_sig = 1.0d0
                    ,amp_fact_init = 0.66d0
                    ,sig_init = 4.0d0
                    ,lb_sig_init = 1.0d0
                    ,ub_sig_init = 12.0d0
                    ,lb_sig = 1.0d0
                    ,ub_sig = 100.0d0
                    ,init_option = 'mean'
                    ,maxiter_init = 15000
                    ,maxiter = 30
                    ,m = 10
                    ,noise = .false.
                    ,regul = .true.
                    ,descent = .true.
                    ,lstd = 1
                    ,ustd = 8
                    ,iprint = -1
                    ,iprint_init = -1
                    ,save_grid = .false.
                    /
                """)
                cmd = addenv(Cmd(`$(joinpath(build, "ROHSA")) $params`; dir=out), "OMP_NUM_THREADS" => string(threads))
                run(pipeline(cmd; stdout=joinpath(out, "stdout.log"), stderr=joinpath(out, "stderr.log")))
                push!(outputs, read(joinpath(out, "out.dat"), String))
            end
            @test !isempty(outputs[1])
            @test outputs[1] == outputs[2]
        end
    end
end
