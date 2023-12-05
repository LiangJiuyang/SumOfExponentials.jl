struct GaussParameter{T}
    sw::Vector{NTuple{2, T}}
end

GaussParameter(Step::Int) = GaussParameter{Float64}([tuple(legendre(Step)[1][i], legendre(Step)[2][i]) for i in 1:Step])

@inline function Gauss_int(integrand::Function, Gaussian::GaussParameter{T}; region::NTuple{2, T} = (T(-1.0), T(1.0))) where {T <: Number}

    a, b = region
    result = sum(i->integrand((b + a) / 2 + (b - a) * i[1] / 2) * (b - a) * i[2] / 2, Gaussian.sw; init= zero(T))
    
    return result
end