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

# Run submit_rohsa.sh with stand-ins for julia (the preflight) and sbatch, and
# return the sbatch arguments it would have submitted with.
function submitted_sbatch_args(config_lines; env=())
    mktempdir() do dir
        bin = joinpath(dir, "bin")
        mkpath(bin)
        write(joinpath(bin, "julia"), "#!/bin/sh\nexit 0\n")
        write(joinpath(bin, "sbatch"), "#!/bin/sh\nprintf '%s\\n' \"\$@\" > \"\$SBATCH_ARGS\"\n")
        chmod(joinpath(bin, "julia"), 0o755)
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
    small = submitted_sbatch_args(["crop_size = 512", "n_gauss = [8]"])
    @test small["--partition"] == "ZEN4"
    @test small["--cpus-per-task"] == "16"
    @test small["--mem"] == "12G"

    large = submitted_sbatch_args(["crop_size = 2048", "n_gauss = [16]"])
    @test large["--cpus-per-task"] == "32"
    @test large["--mem"] == "120G"

    whole_cube = submitted_sbatch_args(["crop_size = 0", "n_gauss = [6]"])
    @test whole_cube["--partition"] == "ZEN4"
    @test whole_cube["--mem"] == "180G"

    override = submitted_sbatch_args(["crop_size = 512", "n_gauss = [8]"]; env=("MEM" => "50G", "CPUS" => "8"))
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
