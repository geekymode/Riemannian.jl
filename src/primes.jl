# ─────────────────────────────────────────────────────────────────────────────
# Chapter 1 — Primes
# ─────────────────────────────────────────────────────────────────────────────

"""
    sieve(n) -> Vector{Int}

All primes `p ≤ n` by the sieve of Eratosthenes.
"""
function sieve(n::Integer)
    n < 2 && return Int[]
    isp = trues(n)
    isp[1] = false
    for p in 2:isqrt(n)
        isp[p] || continue
        for k in p*p:p:n
            isp[k] = false
        end
    end
    return findall(isp)
end

"""
    isprime(n) -> Bool

Primality test: trial division for small `n`, deterministic Miller–Rabin for 64-bit `n`.
"""
function isprime(n::Integer)
    n < 2 && return false
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)
        n == p && return true
        n % p == 0 && return false
    end
    n < 1681 && return true                       # 41² — trial division above was enough
    m = Int128(n)
    d, r = m - 1, 0
    while iseven(d)
        d >>= 1
        r += 1
    end
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)  # deterministic for n < 3.3e24
        x = powermod(Int128(a), d, m)
        (x == 1 || x == m - 1) && continue
        composite = true
        for _ in 1:r-1
            x = mod(x * x, m)
            if x == m - 1
                composite = false
                break
            end
        end
        composite && return false
    end
    return true
end

"""
    primepi_table(n) -> Vector{Int}

Cumulative counts: `primepi_table(n)[k] == π(k)` for `1 ≤ k ≤ n`.
"""
function primepi_table(n::Integer)
    counts = zeros(Int, max(n, 1))
    c = 0
    isp = falses(max(n, 1))
    isp[sieve(n)] .= true
    for k in 1:n
        c += isp[k]
        counts[k] = c
    end
    return counts
end

"""
    primepi(x) -> Int

The prime-counting function π(x) = #{p ≤ x}.
"""
primepi(x::Real) = x < 2 ? 0 : length(sieve(floor(Int, x)))

"""
    nthprime(n) -> Int

The `n`-th prime (`nthprime(1) == 2`), using Rosser's bound p_n < n(log n + log log n).
"""
function nthprime(n::Integer)
    n ≥ 1 || throw(DomainError(n, "n must be ≥ 1"))
    n < 6 && return (2, 3, 5, 7, 11)[n]
    bound = ceil(Int, n * (log(n) + log(log(n))))
    return sieve(bound)[n]
end

"""
    factorize(n) -> Vector{Pair{Int,Int}}

Prime factorisation `[p₁ => e₁, p₂ => e₂, …]` by trial division.
"""
function factorize(n::Integer)
    n ≥ 1 || throw(DomainError(n, "n must be ≥ 1"))
    out = Pair{Int,Int}[]
    m = Int(n)
    for p in Iterators.flatten((2:3, Iterators.countfrom(5, 2)))
        p * p > m && break
        e = 0
        while m % p == 0
            m ÷= p
            e += 1
        end
        e > 0 && push!(out, p => e)
    end
    m > 1 && push!(out, m => 1)
    return out
end

"""
    prime_gaps(n) -> (primes, gaps)

Primes `p ≤ n` (all but the last) together with the gap to the next prime.
"""
function prime_gaps(n::Integer)
    ps = sieve(n)
    return ps[1:end-1], diff(ps)
end

"""
    twin_primes(n) -> Vector{Tuple{Int,Int}}

Twin prime pairs `(p, p+2)` with `p + 2 ≤ n`.
"""
function twin_primes(n::Integer)
    ps = sieve(n)
    return [(ps[i], ps[i+1]) for i in 1:length(ps)-1 if ps[i+1] - ps[i] == 2]
end

# ── classical approximations to π(x) ─────────────────────────────────────────

"""
    li(x)

Logarithmic integral li(x) = PV ∫₀ˣ dt / log t  (= Ei(log x)).  Gauss's guess for π(x).
"""
li(x::Real) = x == 1 ? -Inf : expinti(log(float(x)))

"""
    Li(x)

Offset logarithmic integral Li(x) = li(x) − li(2) = ∫₂ˣ dt / log t.
"""
Li(x::Real) = li(x) - li(2.0)

"""
    riemann_R(x)

Riemann's prime-counting approximation R(x) = Σₙ μ(n)/n · li(x^{1/n}),
evaluated with the rapidly convergent Gram series
R(x) = 1 + Σ_{k≥1} (log x)ᵏ / (k · k! · ζ(k+1)).
"""
function riemann_R(x::Real)
    x ≤ 0 && return 0.0
    L = log(float(x))
    s = 1.0
    term = 1.0                                   # (log x)^k / k!
    for k in 1:500
        term *= L / k
        Δ = term / (k * real(zeta(k + 1.0)))
        s += Δ
        abs(Δ) < 1e-17 * abs(s) && k > L && break
    end
    return s
end
