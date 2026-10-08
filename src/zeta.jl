# ─────────────────────────────────────────────────────────────────────────────
# Chapter 2 — The zeta function: series, Euler product, analytic continuation
# ─────────────────────────────────────────────────────────────────────────────

# ── exact values ─────────────────────────────────────────────────────────────

"""
    zeta_even_exact(n) -> Rational{BigInt}

The rational r with ζ(2n) = r · π^{2n} (Euler, 1735/1740):
ζ(2n) = (−1)^{n+1} B₂ₙ (2π)^{2n} / (2 (2n)!).  `zeta_even_exact(1) == 1//6` (Basel problem).
"""
function zeta_even_exact(n::Integer)
    n ≥ 1 || throw(DomainError(n, "n must be ≥ 1"))
    return (-1)^(n + 1) * bernoulli(2n) * big(2)^(2n - 1) / factorial(big(2n))
end

"""
    zeta_negint_exact(n) -> Rational{BigInt}

ζ(−n) = (−1)ⁿ Bₙ₊₁ / (n+1) for n ≥ 0.  ζ(0) = −1/2, ζ(−1) = −1/12, ζ(−2k) = 0 (trivial zeros).
"""
function zeta_negint_exact(n::Integer)
    n ≥ 0 || throw(DomainError(n, "n must be ≥ 0"))
    return (-1)^n * bernoulli(n + 1) / (n + 1)
end

# ── naive series and product (convergent only for Re s > 1, or Re s > 0 for η) ──

"""
    dirichlet_partial(s, N)

Partial Dirichlet sum Σ_{n=1}^{N} n^{-s}.  Converges to ζ(s) only for Re s > 1.
"""
dirichlet_partial(s::Number, N::Integer) = sum(n -> exp(-s * log(float(n))), 1:N)

"""
    eta_partial(s, N)

Partial alternating sum Σ_{n=1}^{N} (−1)^{n−1} n^{-s}; converges to η(s) for Re s > 0.
"""
eta_partial(s::Number, N::Integer) = sum(n -> (isodd(n) ? 1 : -1) * exp(-s * log(float(n))), 1:N)

"""
    euler_product_partial(s, P)

Partial Euler product Π_{p≤P} (1 − p^{-s})^{-1}.  Equals ζ(s) in the limit for Re s > 1:
the fundamental bridge between primes and ζ.
"""
euler_product_partial(s::Number, P::Integer) = prod(p -> inv(1 - exp(-s * log(float(p)))), sieve(P); init = one(complex(float(s))))

# ── Borwein's algorithm (alternating series acceleration) ────────────────────

"""
    zeta_borwein(s; n = 60)

ζ(s) via P. Borwein's accelerated η-series, ζ(s) = η(s)/(1 − 2^{1−s}).
Accurate for Re s > 0 and modest |Im s| (error grows like e^{π|t|/2}); a nice
illustration that a *convergent* series already continues ζ to 0 < Re s < 1.
"""
function zeta_borwein(s::Number; n::Integer = 60)
    z = complex(float(s))
    T = real(typeof(z))
    d = Vector{T}(undef, n + 1)
    acc = zero(T)
    term = T(1) / n                                 # i = 0 term: (n−1)!/(n! 0!) = 1/n
    for i in 0:n
        if i > 0
            term *= T(4) * (n + i - 1) * (n - i + 1) / ((2i) * (2i - 1))
        end
        acc += term
        d[i+1] = n * acc
    end
    s_ = zero(z)
    for k in 0:n-1
        s_ += (isodd(k) ? -1 : 1) * (d[k+1] - d[n+1]) * exp(-z * log(T(k + 1)))
    end
    η = -s_ / d[n+1]
    return η / (1 - exp((1 - z) * log(T(2))))
end

# ── Euler–Maclaurin: the workhorse ───────────────────────────────────────────

"""
    zeta_em(s; N, M)

ζ(s) by Euler–Maclaurin summation:

ζ(s) = Σ_{n<N} n^{-s} + N^{1−s}/(s−1) + N^{-s}/2 + Σ_{k=1}^{M} B₂ₖ/(2k)! · s(s+1)⋯(s+2k−2) · N^{−s−2k+1} + R.

The right side is analytic for all s ≠ 1, so this formula *is* an analytic continuation.
Defaults pick N, M so the error is near machine precision for the float type of `s`
(works with `BigFloat`). Best for Re s ≥ 0; [`zeta`](@ref) reflects otherwise.
"""
function zeta_em(s::Complex{T}; N::Union{Nothing,Integer} = nothing,
                 M::Union{Nothing,Integer} = nothing) where {T<:AbstractFloat}
    s == 1 && return Complex{T}(Inf, 0)
    bits = precision(T)
    M = something(M, bits ÷ 2 + 5)
    N = something(N, max(10, ceil(Int, (abs(s) + 2M) / π)))
    acc = zero(Complex{T})
    for n in 1:N-1
        acc += exp(-s * log(T(n)))
    end
    logN = log(T(N))
    N⁻ˢ = exp(-s * logN)
    acc += N⁻ˢ * N / (s - 1) + N⁻ˢ / 2
    c = _em_coeffs(T, M)
    rising = s                                      # s(s+1)…(s+2k−2)
    w = N⁻ˢ / N                                     # N^{−s−2k+1} for k = 1
    invN² = inv(T(N)^2)
    for k in 1:M
        acc += c[k] * rising * w
        rising *= (s + 2k - 1) * (s + 2k)
        w *= invN²
    end
    return acc
end
zeta_em(s::Number; kw...) = zeta_em(complex(float(s)); kw...)

# log χ(s) where ζ(s) = χ(s) ζ(1−s), χ(s) = 2ˢ π^{s−1} sin(πs/2) Γ(1−s)
_logchi(s::Complex{T}) where {T} =
    s * log(T(2)) + (s - 1) * log(T(π)) + _logsin(T(π) * s / 2) + loggamma_c(1 - s)

"""
    chi_factor(s)

χ(s) = 2ˢ π^{s−1} sin(πs/2) Γ(1−s), the factor in the functional equation ζ(s) = χ(s) ζ(1−s).
On the critical line |χ(1/2+it)| = 1. The sin(πs/2) factor vanishes at s = −2, −4, …:
these are the *trivial zeros*.
"""
chi_factor(s::Number) = exp(_logchi(complex(float(s))))

"""
    zeta(s)

The Riemann zeta function for any complex `s ≠ 1` (generic float type, incl. `BigFloat`).

* Re s ≥ 0: Euler–Maclaurin ([`zeta_em`](@ref)).
* Re s < 0: functional equation ζ(s) = χ(s) ζ(1−s); exact zero at negative even integers,
  exact Bernoulli values at negative odd integers.

Real input gives a real result.
"""
function zeta(s::Complex{T}) where {T<:AbstractFloat}
    if real(s) < 0
        if imag(s) == 0 && isinteger(real(s))
            n = -Int(real(s))
            return iseven(n) ? zero(Complex{T}) : Complex{T}(T(zeta_negint_exact(n)))
        end
        return exp(_logchi(s)) * zeta_em(1 - s)
    end
    return zeta_em(s)
end
zeta(s::Complex) = zeta(complex(float(real(s)), float(imag(s))))
zeta(s::Real) = s == 1 ? oftype(float(s), Inf) : real(zeta(complex(float(s))))

"""
    dirichlet_eta(s)

Dirichlet eta η(s) = Σ (−1)^{n−1} n^{-s} = (1 − 2^{1−s}) ζ(s), entire.  η(1) = log 2.
"""
function dirichlet_eta(s::Number)
    z = complex(float(s))
    z == 1 && return complex(log(real(oftype(z, 2))))
    v = (1 - exp((1 - z) * log(real(oftype(z, 2))))) * zeta(z)
    return s isa Real ? real(v) : v
end

"""
    completed_zeta(s)

Λ(s) = π^{−s/2} Γ(s/2) ζ(s), satisfying the symmetric functional equation Λ(s) = Λ(1−s).
The Γ(s/2) poles cancel the trivial zeros, so Λ's zeros are exactly the nontrivial ones.
"""
function completed_zeta(s::Number)
    z = complex(float(s))
    T = real(typeof(z))
    return exp(-z / 2 * log(T(π)) + loggamma_c(z / 2)) * zeta(z)
end

"""
    xi(s)

Riemann's ξ(s) = ½ s(s−1) π^{−s/2} Γ(s/2) ζ(s): entire, ξ(s) = ξ(1−s), ξ(s̄) = conj ξ(s),
zeros = nontrivial zeros of ζ. Real on the critical line and on the real axis.
"""
function xi(s::Number)
    z = complex(float(s))
    T = real(typeof(z))
    (z == 0 || z == 1) && return complex(T(0.5))
    return z * (z - 1) / 2 * completed_zeta(z)
end

"""
    Xi(t)

Riemann's Ξ(t) = ξ(1/2 + it), a *real* even function of real t. RH ⇔ all zeros of Ξ are real.
Decays like t^{7/4} e^{−πt/4}, so it underflows in Float64 for t ≳ 900.
"""
Xi(t::Real) = real(xi(complex(oftype(float(t), 0.5), float(t))))

# ── partial sums as paths in the complex plane ───────────────────────────────

"""
    partial_sums(s, N; series = :zeta, corrected = false) -> Vector{Complex}

The path of partial sums Sₙ = Σ_{k≤n} aₖ, n = 1…N, of the Dirichlet series for ζ (`aₖ = k⁻ˢ`)
or η (`aₖ = (−1)^{k−1} k⁻ˢ`). Drawn as a chain of vectors, this is the classic "spiral" picture:

* `:zeta`, raw: for 0 < Re s < 1 the sum diverges, and once k ≳ |t|/2π the path becomes a spiral
  winding *outward* around ζ(s) (Sₙ − ζ(s) ≈ n^{1−s}/(1−s)).
* `:zeta`, `corrected = true`: Sₙ + n^{1−s}/(s−1) − n^{−s}/2 (Euler–Maclaurin) spirals *inward* to ζ(s).
* `:eta`, raw: zig-zags toward η(s) = (1 − 2^{1−s}) ζ(s); `corrected = true` averages consecutive
  sums, which converges much faster. The remaining error still alternates in sign, so plot every
  other point (`[2:2:end]`) for a smooth inward spiral.

At a nontrivial zero the spiral's centre is the origin.
"""
function partial_sums(s::Number, N::Integer; series::Symbol = :zeta, corrected::Bool = false)
    series in (:zeta, :eta) || throw(ArgumentError("series must be :zeta or :eta"))
    z = complex(float(s))
    T = real(typeof(z))
    S = Vector{typeof(z)}(undef, N + 1)
    acc = zero(z)
    for n in 1:N+1
        term = exp(-z * log(T(n)))
        acc += (series === :eta && iseven(n)) ? -term : term
        S[n] = acc
    end
    corrected || return S[1:N]
    if series === :zeta
        return [S[n] + exp((1 - z) * log(T(n))) / (z - 1) - exp(-z * log(T(n))) / 2 for n in 1:N]
    else
        return [(S[n] + S[n+1]) / 2 for n in 1:N]
    end
end
