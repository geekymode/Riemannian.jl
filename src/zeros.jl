# ─────────────────────────────────────────────────────────────────────────────
# Chapter 3 — Zeros: trivial, nontrivial, counting, Gram points
# ─────────────────────────────────────────────────────────────────────────────

"""
    trivial_zeros(n) -> Vector{Int}

The first `n` trivial zeros −2, −4, …, −2n.  They come from the factor sin(πs/2) in the
functional equation (equivalently, from the poles of Γ(s/2) in Λ(s)).
"""
trivial_zeros(n::Integer) = [-2k for k in 1:n]

"""
    KNOWN_ZEROS

Imaginary parts γ of the first ten nontrivial zeros ρ = 1/2 + iγ (A. Odlyzko's tables),
for reference and testing.
"""
const KNOWN_ZEROS = (
    14.134725141734693790, 21.022039638771554993, 25.010857580145688763,
    30.424876125859513210, 32.935061587739189691, 37.586178158825671257,
    40.918719012147495187, 43.327073280914999519, 48.005150881167159727,
    49.773832477672302182,
)

"""
    riemann_siegel_theta(t)

θ(t) = arg Γ(1/4 + it/2) − (t/2) log π  (continuous branch, θ(0) = 0).
Chosen so that e^{iθ(t)} ζ(1/2+it) is real.  θ(t) ≈ (t/2) log(t/2πe) − π/8 + 1/(48t).
"""
function riemann_siegel_theta(t::Real)
    T = float(typeof(t))
    t < 0 && return -riemann_siegel_theta(-t)
    return imag(loggamma_c(complex(T(0.25), t / 2))) - t / 2 * log(T(π))
end

"""
    hardy_Z(t)

Hardy's Z-function Z(t) = e^{iθ(t)} ζ(1/2 + it): real, |Z(t)| = |ζ(1/2+it)|.
Its sign changes are zeros of ζ on the critical line — this is how zeros are *proved*
to lie on the line.
"""
function hardy_Z(t::Real)
    T = float(typeof(t))
    return real(cis(riemann_siegel_theta(T(t))) * zeta(complex(T(0.5), T(t))))
end

"""
    riemann_siegel_Z(t; corrections = true)

Riemann–Siegel formula (found in Riemann's Nachlass by Siegel, 1932):

Z(t) ≈ 2 Σ_{n≤N} n^{−1/2} cos(θ(t) − t log n) + (−1)^{N−1} (t/2π)^{−1/4} C₀(p),

with N = ⌊√(t/2π)⌋, p = frac(√(t/2π)), C₀(p) = cos(2π(p²−p−1/16))/cos(2πp).
Only O(√t) terms instead of O(t): this is what makes large-height computations feasible.
Error is O(t^{−3/4}).
"""
function riemann_siegel_Z(t::Real; corrections::Bool = true)
    t = float(t)
    a = sqrt(t / 2π)
    N = floor(Int, a)
    θ = riemann_siegel_theta(t)
    s = 0.0
    for n in 1:N
        s += cos(θ - t * log(n)) / sqrt(n)
    end
    Z = 2s
    if corrections
        p = a - N
        Z += (-1)^(N - 1) * (t / 2π)^(-1 / 4) * _rs_C0(p)
    end
    return Z
end

function _rs_C0(p)
    c = cos(2π * p)
    if abs(c) < 1e-6                     # removable singularities at p = 1/4, 3/4
        h = 1e-4
        return (_rs_C0(p - h) + _rs_C0(p + h)) / 2
    end
    return cos(2π * (p^2 - p - 1 / 16)) / c
end

# ── Gram points ──────────────────────────────────────────────────────────────

"""
    gram_point(n)

The n-th Gram point gₙ: θ(gₙ) = nπ (g₀ ≈ 17.8456).
"""
function gram_point(n::Integer)
    t = 2π * exp(1 + _lambertw((8n + 1) / (8ℯ)))       # asymptotic inverse of θ
    for _ in 1:50
        Δ = (riemann_siegel_theta(t) - n * π) / (log(t / 2π) / 2)
        t -= Δ
        abs(Δ) < 1e-13 * t && break
    end
    return t
end

"""
    gram_points(n1, n2) -> Vector{Float64}
"""
gram_points(n1::Integer, n2::Integer) = [gram_point(n) for n in n1:n2]

"""
    gram_law_violations(nmax) -> Vector{Int}

Indices n ≤ nmax where Gram's law (−1)ⁿ Z(gₙ) > 0 fails. The first is n = 126.
"""
gram_law_violations(nmax::Integer) = [n for n in 0:nmax if (-1)^n * hardy_Z(gram_point(n)) ≤ 0]

# ── counting zeros ───────────────────────────────────────────────────────────

"""
    argzeta_S(t)

S(t) = π⁻¹ arg ζ(1/2 + it), with arg defined by continuous variation along the horizontal
line from +∞ + it (arg = 0 there) to 1/2 + it.  Small and wildly oscillating: S(t) = O(log t).
"""
function argzeta_S(t::Real; σ0 = 3.0)
    σ = float(σ0)
    z = zeta(complex(σ, t))
    arg_ = angle(z)                                  # principal value is right at σ0 = 3
    h = 0.05
    while σ > 0.5
        hσ = min(h, σ - 0.5)
        znew = zeta(complex(σ - hσ, t))
        Δ = angle(znew / z)
        if abs(Δ) > π / 4 && hσ > 1e-6
            h = hσ / 2
            continue
        end
        arg_ += Δ
        σ -= hσ
        z = znew
        h = min(2h, 0.05)
    end
    return arg_ / π
end

"""
    riemann_von_mangoldt(T)

Smooth part of the zero count: θ(T)/π + 1 ≈ (T/2π) log(T/2πe) + 7/8.
"""
riemann_von_mangoldt(T::Real) = riemann_siegel_theta(T) / π + 1

"""
    zero_count(T) -> Int

N(T) = #{ρ : 0 < Im ρ ≤ T}, computed rigorously-in-spirit via N(T) = θ(T)/π + 1 + S(T)
(argument principle).  Comparing with the number of sign changes of Z on (0, T]
verifies that *all* zeros up to height T lie on the critical line.
"""
zero_count(T::Real) = round(Int, riemann_von_mangoldt(T) + argzeta_S(T))

# ── locating zeros ───────────────────────────────────────────────────────────

# Illinois (modified regula falsi) root refinement on a bracketing interval.
function _refine(f, a, b, fa, fb; tol = 1e-13)
    side = 0
    c = a
    for _ in 1:100
        c = (a * fb - b * fa) / (fb - fa)
        fc = f(c)
        (abs(b - a) < tol * max(1, abs(c)) || fc == 0) && return c
        if fc * fb > 0
            b, fb = c, fc
            side == -1 && (fa /= 2)
            side = -1
        else
            a, fa = c, fc
            side == +1 && (fb /= 2)
            side = +1
        end
    end
    return c
end

"""
    zeros_between(t1, t2; refine = 4, verify = true) -> Vector{Float64}

Heights γ ∈ (t1, t2] of zeros of ζ on the critical line, found as sign changes of `hardy_Z`
on a grid of ~`refine` points per average gap 2π/log(t/2π), then refined to ~1e-13.

With `verify = true` the count is compared with N(t2) − N(t1) from [`zero_count`](@ref);
on mismatch (e.g. a close "Lehmer pair") the grid is refined; a warning is raised if it
still disagrees.
"""
function zeros_between(t1::Real, t2::Real; refine::Integer = 4, verify::Bool = true, maxlevel = 5)
    expected = verify ? zero_count(t2) - zero_count(t1) : -1
    found = Float64[]
    for level in 0:maxlevel
        found = _scan_zeros(float(t1), float(t2), refine * 2^level)
        (!verify || length(found) == expected) && return found
    end
    @warn "zeros_between($t1, $t2): found $(length(found)) sign changes, N(T) predicts $expected"
    return found
end

function _scan_zeros(t1, t2, perGap)
    out = Float64[]
    t = t1
    Zt = hardy_Z(t)
    while t < t2
        gap = 2π / log(max(t, 20.0) / 2π)
        tn = min(t + gap / perGap, t2)
        Zn = hardy_Z(tn)
        if Zt * Zn < 0
            push!(out, _refine(hardy_Z, t, tn, Zt, Zn))
        elseif Zn == 0
            push!(out, tn)
        end
        t, Zt = tn, Zn
    end
    return out
end

const _ZERO_CACHE = Float64[]

"""
    nontrivial_zeros(n) -> Vector{Float64}

Heights γ₁ < γ₂ < … < γₙ of the first `n` nontrivial zeros ρ = 1/2 + iγ (computed, cached).
The first few thousand take seconds.
"""
function nontrivial_zeros(n::Integer)
    while length(_ZERO_CACHE) < n
        lo = isempty(_ZERO_CACHE) ? 10.0 : last(_ZERO_CACHE) + 1e-4
        # blocks of ~100 zeros (or fewer if fewer are needed), each verified against N(T),
        # so a close pair only forces a finer grid locally
        need = min(n - length(_ZERO_CACHE), 100)
        hi = lo + max(20.0, 1.1 * need * 2π / log(max(lo, 20.0) / 2π))
        append!(_ZERO_CACHE, zeros_between(lo, hi))
    end
    return _ZERO_CACHE[1:n]
end

"""
    check_zeros(T) -> NamedTuple

Turing-style verification up to height T: sign changes of Z vs. N(T) from the argument principle.
"""
function check_zeros(T::Real)
    found = zeros_between(10.0, T; verify = false, refine = 8)
    N = zero_count(T)
    return (T = T, sign_changes = length(found), N_T = N, all_on_line = length(found) == N)
end

"""
    close_pairs(γs; k = 10) -> Vector{NamedTuple}

The `k` closest consecutive pairs of zeros measured in units of mean spacing — the
"Lehmer pairs" (near-misses of RH; e.g. near t ≈ 7005.06) relevant to de Bruijn–Newman.
"""
function close_pairs(γs::AbstractVector; k::Integer = 10)
    δ = normalized_spacings(γs)
    idx = partialsortperm(δ, 1:min(k, length(δ)))
    return [(γ1 = γs[i], γ2 = γs[i+1], spacing = δ[i]) for i in idx]
end

# ── how wide is the band where zeros can be? ─────────────────────────────────

"""
    RH_VERIFIED_HEIGHT

All nontrivial zeros with 0 < Im ρ ≤ 3·10¹² lie on the critical line (Platt & Trudgian, 2021).
Below this height the "band" of possible zeros has width 0.
"""
const RH_VERIFIED_HEIGHT = 3.0e12

# Width δ of the proven zero-free sliver σ > 1 − δ, as a function of L = log t (avoids overflow).
function _zero_free_width(L::Real, method::Symbol)
    if method === :classical
        return 1 / (5.573412 * L)
    elseif method === :vinogradov_korobov
        return 1 / (57.54 * L^(2 / 3) * log(L)^(1 / 3))
    else
        throw(ArgumentError("method must be :classical or :vinogradov_korobov"))
    end
end

"""
    zero_free_boundary(t; method = :classical) -> σ₀

Proven zero-free region: ζ(σ + it) ≠ 0 for σ ≥ σ₀(t). By the functional equation also for σ ≤ 1 − σ₀(t).

* `:classical` — de la Vallée Poussin shape, explicit constant of Mossinghoff & Trudgian (2015):
  σ₀ = 1 − 1/(5.573412 log t), valid for t ≥ 2.
* `:vinogradov_korobov` — explicit Vinogradov–Korobov region of Ford (2002):
  σ₀ = 1 − 1/(57.54 (log t)^{2/3} (log log t)^{1/3}), valid for t ≥ 3; better only for log t ≳ 10⁴.

The sliver width 1 − σ₀ → 0 as t → ∞: what we can *prove* about the band of possible zeros
approaches the full critical strip of width 1. RH says the band is the single line Re s = 1/2.
"""
zero_free_boundary(t::Real; method::Symbol = :classical) = 1 - _zero_free_width(log(t), method)
