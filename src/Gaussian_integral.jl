struct GaussParameter{T}
    sw::Vector{NTuple{2, T}}
end

function GaussParameter(Step::Int)  
    legendre_sw = legendre(Step)
    return GaussParameter{Float64}([tuple(legendre_sw[1][i], legendre_sw[2][i]) for i in 1:Step])
end

function GaussParameter(T::DataType, Step::Int) 
    legendre_sw = legendre(T, Step)
    return GaussParameter{T}([tuple(legendre_sw[1][i], legendre_sw[2][i]) for i in 1:Step])
end

@inline function Gauss_int(integrand::Function, Gaussian::GaussParameter{TG}; region::NTuple{2, T} = (T(-1.0), T(1.0))) where {T <: Number, TG <: Number}

    a, b = TG.(region)
    result = sum(i->integrand((b + a) / TG(2) + (b - a) * i[1] / TG(2)) * (b - a) * i[2] / TG(2), Gaussian.sw; init= zero(TG))
    
    return result
end