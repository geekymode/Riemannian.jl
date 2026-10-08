# ─────────────────────────────────────────────────────────────────────────────
# Chapter 7 — The de Bruijn–Newman constant Λ (heat flow on Ξ)
#
#   H_t(z) = ∫₀^∞ e^{t u²} Φ(u) cos(z u) du,   H₀(z) = Ξ(z/2) / 8.
#
# There is a constant Λ with: H_t has only real zeros ⇔ t ≥ Λ.  RH ⇔ Λ ≤ 0.
# Newman (1976) conjectured Λ ≥ 0 ("RH, if true, is barely true"),
# proved by Rodgers & Tao (2018/2020).  Polymath15 (2019): Λ ≤ 0.22;
# Platt & Trudgian (2021): Λ ≤ 0.2.  So 0 ≤ Λ ≤ 0.2.
# ─────────────────────────────────────────────────────────────────────────────

"""
    debruijn_Phi(u)

Φ(u) = Σ_{n≥1} (2π²n⁴e^{9u} − 3πn²e^{5u}) exp(−πn²e^{4u}), the (super-exponentially decaying,
even) Fourier kernel of Ξ.
"""
function debruijn_Phi(u::T) where {T<:AbstractFloat}
    πT = T(π)
    e4 = exp(4u)
    s = zero(T)
    for n in 1:20
        g = πT * n^2 * e4
        g > precision(T) * 0.7 + 50 && break
        s += (2πT^2 * n^4 * exp(9u) - 3πT * n^2 * exp(5u)) * exp(-g)
    end
    return s
end
debruijn_Phi(u::Real) = debruijn_Phi(float(u))

const _PHI_GRID = Dict{Tuple{Int,Float64},Vector{BigFloat}}()

# Φ(kh), k = 0, 1, … until negligible; independent of (t, z), so cached per (prec, h).
function _phi_grid(prec::Integer, h::Real)
    get!(_PHI_GRID, (Int(prec), Float64(h))) do
        setprecision(BigFloat, prec) do
            H = BigFloat(h)
            out = BigFloat[]
            k = 0
            while true
                φ = debruijn_Phi(k * H)
                push!(out, φ)
                (k * H > 0.5 && abs(φ) < big(2.0)^(-prec - 20)) && break
                k += 1
            end
            out
        end
    end
end

"""
    debruijn_H(t, z; prec = 128, h = 1/200)

H_t(z) by the trapezoidal rule (exponentially accurate for this analytic, rapidly
decaying integrand), in `BigFloat` with `prec` bits because H_t(z) ~ e^{−πz/8} is
tiny compared with the integrand.  Returns a `Float64`.
"""
function debruijn_H(t::Real, z::Real; prec::Integer = 128, h::Real = 1 / 200)
    Φ = _phi_grid(prec, h)
    setprecision(BigFloat, prec) do
        T, Z, H = BigFloat(t), BigFloat(z), BigFloat(h)
        s = Φ[1] / 2
        for k in 1:length(Φ)-1
            u = k * H
            s += exp(T * u^2) * Φ[k+1] * cos(Z * u)
        end
        Float64(s * H)
    end
end

"""
    heat_flow_zeros(t; zmax = 100, step = 0.25, prec = 128) -> Vector{Float64}

Real zeros of H_t in (0, zmax].  For t = 0 these are 2γₙ (twice the zeta-zero heights).
As t increases zeros repel and spread; as t decreases they attract and eventually collide
and leave the real axis — Λ is the last time any collision happens.
"""
function heat_flow_zeros(t::Real; zmax::Real = 100, step::Real = 0.25, prec::Integer = 128)
    f(x) = debruijn_H(t, x; prec)
    out = Float64[]
    x = step
    fx = f(x)
    while x < zmax
        xn = min(x + step, zmax)
        fn = f(xn)
        if fx * fn < 0
            push!(out, _refine(f, x, xn, fx, fn; tol = 1e-10))
        end
        x, fx = xn, fn
    end
    return out
end
