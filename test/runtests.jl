using VGT
using RohsaLMC
using Test

@testset "ROHSA-LMC" begin
    include("rohsa_auto_config.jl")
    include("rohsa_products.jl")
    include("rohsa_memory.jl")
    include("rohsa_threads.jl")
    include("rohsa_noise.jl")
    include("rohsa_layers.jl")
    include("derotate.jl")
    include("rohsa_hpc.jl")
    include("script_smoke.jl")
end
