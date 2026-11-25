using FITSIO

module RohsaIrgetasPreflightTest
    include(joinpath(@__DIR__, "..", "hpc", "irgetas", "rohsa_preflight.jl"))
end

@testset "Irgetas ROHSA resource classification" begin
    @test RohsaIrgetasPreflightTest.irgetas_partition(100.0) == :ZEN4
    @test RohsaIrgetasPreflightTest.irgetas_partition(400.0) == :H100
    @test RohsaIrgetasPreflightTest.irgetas_partition(800.0) == :none

    slurm = read(joinpath(@__DIR__, "..", "hpc", "irgetas", "rohsa.slurm"), String)
    @test occursin("#SBATCH --nodes=1", slurm)
    @test occursin("#SBATCH --ntasks=1", slurm)
    @test occursin("srun --ntasks=1", slurm)
    @test !occursin("--gres=gpu", slurm)
end

# Run submit_rohsa.sh with a stand-in for sbatch, and return the sbatch
# arguments it would have submitted with. By default julia is a stand-in that
# reports max_n as the preflight's largest candidate N; with real_julia the
# real preflight runs on the config (which must be complete).
function submitted_sbatch_args(config_lines; env=(), max_n=6, real_julia=false)
    mktempdir() do dir
        bin = joinpath(dir, "bin")
        mkpath(bin)
        real_julia || write(joinpath(bin, "julia"), "#!/bin/sh\necho 'ROHSA largest candidate N: $max_n'\n")
        write(joinpath(bin, "sbatch"), "#!/bin/sh\nprintf '%s\\n' \"\$@\" > \"\$SBATCH_ARGS\"\n")
        real_julia || chmod(joinpath(bin, "julia"), 0o755)
        chmod(joinpath(bin, "sbatch"), 0o755)
        config = joinpath(dir, "run.toml")
        write(config, join(config_lines, "\n") * "\n")
        args_file = joinpath(dir, "sbatch_args")
        script = joinpath(@__DIR__, "..", "hpc", "irgetas", "submit_rohsa.sh")
        path = bin * ":" * ENV["PATH"]
        cmd = setenv(`bash $script $config`, Dict("PATH" => path, "SBATCH_ARGS" => args_file, env...); dir)
        run(pipeline(cmd; stdout=devnull))
        return Dict(split(arg, "=", limit=2)[1] => get(split(arg, "=", limit=2), 2, "")
                    for arg in readlines(args_file) if startswith(arg, "--"))
    end
end

@testset "submit_rohsa.sh sizes the Slurm job from the config" begin
    small = submitted_sbatch_args(["crop_size = 512", "n_gauss = [8]"]; max_n=8)
    @test small["--partition"] == "ZEN4"
    @test small["--cpus-per-task"] == "16"
    @test small["--mem"] == "12G"

    large = submitted_sbatch_args(["crop_size = 2048", "n_gauss = [16]"]; max_n=16)
    @test large["--cpus-per-task"] == "32"
    @test large["--mem"] == "120G"

    whole_cube = submitted_sbatch_args(["crop_size = 0", "n_gauss = [6]"])
    @test whole_cube["--partition"] == "ZEN4"
    @test whole_cube["--mem"] == "180G"

    for (n, mem) in ((6, "100G"), (7, "140G"), (8, "140G"), (9, "120G"))
        @test submitted_sbatch_args(["crop_size = 2048"]; max_n=n)["--mem"] == mem
    end

    override = submitted_sbatch_args(["crop_size = 512", "n_gauss = [8]"]; max_n=8, env=("MEM" => "50G", "CPUS" => "8"))
    @test override["--mem"] == "50G"
    @test override["--cpus-per-task"] == "8"
end

# inspect_config on a tiny cube: the estimate must cover the largest candidate N,
# which the [auto_config] n_gauss axis can raise above [rohsa] n_gauss.
function preflight_estimate(auto_axis)
    mktempdir() do dir
        input = joinpath(dir, "cube.fits")
        f = FITS(input, "w")
        try
            write(f, zeros(Float32, 16, 16, 8))
        finally
            close(f)
        end
        toml = read(joinpath(@__DIR__, "..", "configs", "rohsa", "full_cube.toml"), String)
        toml = replace(toml, r"^input = .*$"m => "input = \"$input\"")
        isempty(auto_axis) || (toml = replace(toml, "[auto_config]\nenabled = false\n" =>
                                              "[auto_config]\nenabled = true\nn_gauss = $auto_axis\n"))
        config = joinpath(dir, "run.toml")
        write(config, toml)
        return RohsaIrgetasPreflightTest.inspect_config(config)
    end
end

@testset "ROHSA preflight sizes from the largest candidate N" begin
    est(n) = RohsaIrgetasPreflightTest.rohsa_memory_estimate(16, 16, 8, n, 10)
    plain = preflight_estimate("")
    @test plain.estimate == est(6)
    swept = preflight_estimate("[4, 6, 8]")
    @test swept.estimate == est(8)
    @test swept.estimate.peak_gib > plain.estimate.peak_gib
end

# End to end: the real preflight resolves the [auto_config] n_gauss axis, so a
# sweep to N = 8 is requested at the 7-8 bracket (140G), not at the [rohsa]
# value N = 6 (100G).
@testset "submit_rohsa.sh requests memory for the largest candidate N" begin
    got = mktempdir() do dir
        input = joinpath(dir, "cube.fits")
        f = FITS(input, "w")
        try
            write(f, zeros(Float32, 16, 16, 8))
        finally
            close(f)
        end
        toml = read(joinpath(@__DIR__, "..", "configs", "rohsa", "full_cube.toml"), String)
        toml = replace(toml, r"^input = .*$"m => "input = \"$input\"", r"^crop_size = .*$"m => "crop_size = 2048")
        plain = submitted_sbatch_args(split(toml, "\n"); real_julia=true)
        sweep_toml = replace(toml, "[auto_config]\nenabled = false\n" => "[auto_config]\nenabled = true\nn_gauss = [4, 6, 8]\n")
        sweep = submitted_sbatch_args(split(sweep_toml, "\n"); real_julia=true)
        (plain=plain["--mem"], sweep=sweep["--mem"])
    end
    @test got.plain == "100G"
    @test got.sweep == "140G"
end
