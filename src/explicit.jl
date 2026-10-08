# ─────────────────────────────────────────────────────────────────────────────
# Chapter 4 — Explicit formulas: primes from zeros
# ─────────────────────────────────────────────────────────────────────────────

# Exponential integral Ei for complex argument off the real axis.
_Ei(z::Complex) = -expint(-z) + (imag(z) > 0 ? im * π : -im * π)

"""
    psi_explicit(x, γs)

von Mangoldt's explicit formula truncated to the zeros 1/2 ± iγ for γ ∈ `γs`:

ψ₀(x) = x − Σ_ρ x^ρ/ρ − log 2π − ½ log(1 − x⁻²).

Each conjugate pair contributes a wave 2 Re(x^ρ/ρ) of amplitude ~ √x/|ρ|:
RH says every wave has amplitude √x, which is why the error in the PNT would be O(√x log² x).
"""
function psi_explicit(x::Real, γs::AbstractVector)
    x ≤ 1 && return 0.0
    L = log(x)
    s = 0.0
    for γ in γs
        ρ = complex(0.5, γ)
        s += 2 * real(exp(ρ * L) / ρ)
    end
    return x - s - log(2π) - log1p(-1 / x^2) / 2
end

"""
    riemann_J(x)

Riemann's prime-power counting function J(x) = Σ_{pᵏ≤x} 1/k = Σₙ π(x^{1/n})/n.
"""
function riemann_J(x::Real)
    x < 2 && return 0.0
    s = 0.0
    n = 1
    while x^(1 / n) ≥ 2
        s += primepi(x^(1 / n)) / n
        n += 1
    end
    return s
end

# ∫ₓ^∞ dt / (t (t²−1) log t) = ∫_{log x}^∞ du / ((e^{2u} − 1) u), Simpson's rule.
function _J_tail(x::Real)
    a = log(x)
    b = a + 30
    n = 600
    h = (b - a) / n
    f(u) = 1 / (expm1(2u) * u)
    s = f(a) + f(b)
    for k in 1:n-1
        s += (isodd(k) ? 4 : 2) * f(a + k * h)
    end
    return s * h / 3
end

"""
    riemann_J_explicit(x, γs)

Riemann's 1859 explicit formula, truncated to the given zeros:

J(x) = li(x) − Σ_ρ li(x^ρ) − log 2 + ∫ₓ^∞ dt / (t(t²−1) log t),

with li(x^ρ) interpreted as Ei(ρ log x).
"""
function riemann_J_explicit(x::Real, γs::AbstractVector)
    x ≤ 1 && return 0.0
    L = log(x)
    s = 0.0
    for γ in γs
        s += 2 * real(_Ei(complex(0.5, γ) * L))
    end
    return li(x) - s - log(2.0) + _J_tail(x)
end

"""
    primepi_explicit(x, γs)

π(x) reconstructed from zeros by Möbius inversion: π(x) = Σₙ μ(n)/n · J(x^{1/n}).
With more zeros the smooth curve sharpens into the prime staircase.
"""
function primepi_explicit(x::Real, γs::AbstractVector)
    x < 2 && return 0.0
    s = 0.0
    n = 1
    while x^(1 / n) ≥ 2
        μ = mobius(n)
        μ != 0 && (s += μ / n * riemann_J_explicit(x^(1 / n), γs))
        n += 1
    end
    return s
end
