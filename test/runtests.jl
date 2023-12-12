using SumOfExponentials
using Test
using SpecialFunctions

@testset "SumOfExponentials.jl" begin
    include("VP.jl")
    include("VPMR.jl")
    include("fast_exp_transform.jl")
end
