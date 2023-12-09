module VPMR

using LinearAlgebra, SpecialFunctions, GaussQuadrature, SpecialFunctions, MatrixEquations, GenericLinearAlgebra, GenericSchur

export GaussParameter, Gauss_int
export VP_cal, MR_cal, VPMR_cal
export soe, soe_error

include("Gaussian_integral.jl")
include("VP.jl")
include("MR.jl")
include("soe.jl")

function VPMR_cal(f::Function, nc::T, n::Int, gausspara::GaussParameter{TG}, p::Int; T1::DataType = ComplexF64, T2::DataType = Float64, digit::Int = 1024) where{T, TG}

    @assert iszero(digit % 256)

    s, w = setprecision(digit) do
        VP(f, nc, n, gausspara)
    end

    @show soe_error(f, s, w)

    smr, wmr, σ = setprecision(digit) do 
        MR(s, w, p)
    end

    @show soe_error(f, smr, wmr)

    return T1.(smr), T1.(wmr), T2(σ)
end

end
