# ─────────────────────────────────────────────────────────────────────────────
# Arithmetic functions that encode the primes
# ─────────────────────────────────────────────────────────────────────────────

"""
    mobius(n)

Möbius function μ(n): 0 if n has a squared prime factor, else (−1)^(number of primes).
1/ζ(s) = Σ μ(n) n^{-s}.
"""
function mobius(n::Integer)
    f = factorize(n)
    any(p -> p.second > 1, f) && return 0
    return isodd(length(f)) ? -1 : 1
end

"""
    mobius_table(N) -> Vector{Int}

μ(1), …, μ(N) by a linear sieve.
"""
function mobius_table(N::Integer)
    μ = ones(Int, N)
    composite = falses(N)
    primes = Int[]
    for i in 2:N
        if !composite[i]
            push!(primes, i)
            μ[i] = -1
        end
        for p in primes
            k = i * p
            k > N && break
            composite[k] = true
            if i % p == 0
                μ[k] = 0
                break
            end
            μ[k] = -μ[i]
        end
    end
    return μ
end

"""
    liouville(n)

Liouville function λ(n) = (−1)^Ω(n), Ω = number of prime factors with multiplicity.
ζ(2s)/ζ(s) = Σ λ(n) n^{-s}.
"""
liouville(n::Integer) = iseven(sum(last, factorize(n); init = 0)) ? 1 : -1

"""
    mertens_table(N) -> Vector{Int}

Mertens function M(x) = Σ_{n≤x} μ(n) for x = 1:N.
RH ⇔ M(x) = O(x^{1/2+ε}).
"""
mertens_table(N::Integer) = cumsum(mobius_table(N))

"""
    mertens(x)

Mertens function M(x) = Σ_{n≤x} μ(n).
"""
mertens(x::Real) = x < 1 ? 0 : last(mertens_table(floor(Int, x)))

"""
    vonmangoldt(n)

von Mangoldt Λ(n) = log p if n = pᵏ, else 0.   −ζ'/ζ(s) = Σ Λ(n) n^{-s}.
"""
function vonmangoldt(n::Integer)
    n < 2 && return 0.0
    f = factorize(n)
    return length(f) == 1 ? log(float(f[1].first)) : 0.0
end

"""
    chebyshev_theta(x)

ϑ(x) = Σ_{p≤x} log p.  The prime number theorem is ϑ(x) ~ x.
"""
chebyshev_theta(x::Real) = sum(log ∘ float, sieve(floor(Int, max(x, 0))); init = 0.0)

"""
    chebyshev_psi(x)

ψ(x) = Σ_{n≤x} Λ(n) = Σ_{pᵏ≤x} log p.  The explicit formula expresses ψ through the zeros of ζ.
"""
function chebyshev_psi(x::Real)
    x < 2 && return 0.0
    X = floor(Int, x)
    s = 0.0
    for p in sieve(X)
        q = p
        while q ≤ X
            s += log(float(p))
            q > X ÷ p && break
            q *= p
        end
    end
    return s
end

"""
    totient(n)

Euler's φ(n) = n Π_{p|n} (1 − 1/p).
"""
function totient(n::Integer)
    r = n
    for (p, _) in factorize(n)
        r = r ÷ p * (p - 1)
    end
    return r
end

"""
    divisor_sigma(n, k = 1)

σₖ(n) = Σ_{d|n} dᵏ.  σ = σ₁ appears in Robin's and Lagarias' RH criteria.
"""
function divisor_sigma(n::Integer, k::Integer = 1)
    s = one(widen(n))
    for (p, e) in factorize(n)
        pk = widen(p)^k
        s *= k == 0 ? (e + 1) : (pk^(e + 1) - 1) ÷ (pk - 1)
    end
    return s
end

"""
    sigma_table(N) -> Vector{Int}

σ(1), …, σ(N) by a divisor sieve.
"""
function sigma_table(N::Integer)
    σ = zeros(Int, N)
    for d in 1:N, m in d:d:N
        σ[m] += d
    end
    return σ
end
