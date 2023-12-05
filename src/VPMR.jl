module VPMR

using LinearAlgebra, SpecialFunctions, GaussQuadrature

export GaussParameter, Gauss_int, VP_kernel, VP_k

include("Gaussian_integral.jl")
include("VP_kernel.jl")


end
