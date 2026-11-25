@testset "vendored ROHSA is the upstream base plus the patch series" begin
    script = joinpath(@__DIR__, "..", "scripts", "check_rohsa_patches.sh")
    code = run(ignorestatus(`bash $script`)).exitcode
    if code == 3
        @test_skip code == 0   # shallow clone: the base commit is not available
    else
        @test code == 0
    end
end
