function VP_kernel(f::Function, nc::T, j::Int) where{T}
    VP_K = r -> f(- nc * log((T(1.0) + cos(T(r))) / T(2))) * cos(T(j) * T(r))
    return VP_K
end

function VP_k(f::Function, nc::TN, n::Int, j::Int, gausspara::GaussParameter{T}) where{TN, T}
    nc = T(nc)

    K = VP_kernel(f, nc, j)

    ∫K = Gauss_int(K, gausspara; region = (T(0.0), T(π)))

    kj = max(T(2)/T(π), (T(4) * n - T(2) * j) / (n * T(π))) * ∫K

    return kj
end

function sign_bit(l::BigInt)
    return l % 2 == 0 ? big(1) : big(-1)
end

function VP_cn(n::BigInt)
    cn = Vector{Vector{BigFloat}}()

    for j in 0 : 2*n-1
        cnj = Vector{BigFloat}()

        for l in 1 : n
            cnjl = sign_bit(n + l - j) * (big(1.0) - l / n) * (n + l) / (n + l + j) * binomial(n + l + j, n + l - j)
            push!(cnj, cnjl)
        end

        push!(cn, cnj)
    end

    return cn
end

function get_cn(cn::Vector{Vector{BigFloat}}, j::T, l::T) where{T<:Integer}
    return cn[j + 1][l]
end

function VP(f::Function, nc::T, n::Int, N::Int) where{T}
    gausspara = GaussParameter(BigFloat, N)

    k0 = VP_k(f, nc, n, 0, gausspara)
    k = big.([VP_k(f, nc, n, j, gausspara) for j in 1:2 * n - 1])

    s_array = [big(j) / big(nc) for j in 1:2*n - 1]
    w_array = Vector{BigFloat}()
    
    n = BigInt(n)

    cn = VP_cn(n)

    # j = 0
    w0 = 2 * k0 + sum([sign_bit(l) * n / (2 * n - l) * k[l] for l in 1:n]) + sum([sign_bit(n + l) * (n - l) / (2 * n - l) * k[n + l] for l in 1:n - 1])

    # push!(w_array, w0)

    #1 ≤ j ≤ n
    for j in 1:n
        wj = sum([sign_bit(l - j) * n * l / ((l + j) * (2 * n - l)) * binomial(l + j, l - j) * k[l] for l in j:n]) 
        wj += sum([get_cn(cn, j, l) * n / (2 * n - l) * k[n + l] for l in 1:n-1])

        wj *= (big(4.0))^j

        push!(w_array, wj)
    end

    #n + 1 ≤ j ≤ 2n - 1
    for j in n + 1:2 * n - 1
        wj = sum([n * get_cn(cn, j, l) * k[n + l] / (2 * n - l) for l in j-n:n-1])
        wj *= (big(4.0))^j
        push!(w_array, wj)
    end

    # return [(s_array[i], w_array[i]) for i in 1:length(s_array)]
    return s_array, w_array
end

function VP_cal(f::Function, nc::T, n::Int, N::Int; digit::Int = 1024) where{T}
    @assert iszero(digit % 256)
    s_array, w_array = setprecision(digit) do 
        VP(f, nc, n, N)
    end
    return s_array, w_array
end