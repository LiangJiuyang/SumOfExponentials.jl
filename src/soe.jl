struct SoePara{T} 
    sw::Vector{Tuple{T, T}}
end

function soe(x::T1, sw::SoePara{T2}) where{T1<:Real, T2}
    x = abs(x)
    sum = zero(ComplexF64)
    for (s, w) in sw.sw
        sum += w * exp(- s * x)
    end
    return T1(real(sum))
end

function soe(x::T, s::Vector{T1}, w::Vector{T2}) where{T, T1, T2}
    return sum([w[i] * exp(-s[i] * x) for i in 1:length(s)])
end

function soe_error(f::Function, s::Vector{T2}, w::Vector{T2}; x::Vector{T1} = big.([0.0:0.01:10.0...])) where{T1<:Real, T2}
    error = [abs(soe(x[i], s, w) - f(x[i])) for i in 1:size(x, 1)]
    return maximum(error)
end