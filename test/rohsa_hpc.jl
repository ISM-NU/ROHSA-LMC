module RohsaIrgetasPreflightTest
    include(joinpath(@__DIR__, "..", "rohsa", "hpc", "irgetas", "preflight.jl"))
end

@testset "Irgetas ROHSA resource classification" begin
    @test RohsaIrgetasPreflightTest.irgetas_partition(100.0) == :ZEN4
    @test RohsaIrgetasPreflightTest.irgetas_partition(400.0) == :H100
    @test RohsaIrgetasPreflightTest.irgetas_partition(800.0) == :none

    slurm = read(joinpath(@__DIR__, "..", "rohsa", "hpc", "irgetas", "run.slurm"), String)
    @test occursin("#SBATCH --nodes=1", slurm)
    @test occursin("#SBATCH --ntasks=1", slurm)
    @test occursin("srun --ntasks=1", slurm)
    @test !occursin("--gres=gpu", slurm)

    submit = read(joinpath(@__DIR__, "..", "rohsa", "hpc", "irgetas", "submit.sh"), String)
    @test occursin("partition=\"ZEN4\"", submit)
    @test occursin("memory=\"180G\"", submit)
end
