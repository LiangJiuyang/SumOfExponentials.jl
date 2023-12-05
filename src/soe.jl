function soe(x::T, s::Vector{T1}, w::Vector{T2}) where{T, T1, T2}
    return sum([w[i] * exp(-s[i] * x) for i in 1:length(s)])
end