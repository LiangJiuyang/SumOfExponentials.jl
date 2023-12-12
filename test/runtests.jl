using SOEVPMR
using Test
using SpecialFunctions

@testset "SOEVPMR.jl" begin
    include("VP.jl")
    include("VPMR.jl")
end
