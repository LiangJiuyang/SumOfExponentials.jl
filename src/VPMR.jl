module VPMR

using LinearAlgebra, SpecialFunctions, GaussQuadrature

export GaussParameter, Gauss_int, VP_sw, soe

include("Gaussian_integral.jl")
include("VP_kernel.jl")
include("soe.jl")

end
