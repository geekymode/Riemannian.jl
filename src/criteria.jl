# ─────────────────────────────────────────────────────────────────────────────
# Chapter 6 — Statements equivalent to RH, and bounds towards it
# ─────────────────────────────────────────────────────────────────────────────

"""
    robin_ratio(n)

σ(n) / (n log log n).  Robin (1984): RH ⇔ robin_ratio(n) < e^γ ≈ 1.78107 for every n > 5040.
"""
robin_ratio(n::Integer) = divisor_sigma(n) / (n * log(log(n)))

"""
    robin_violations(N) -> Vector{Int}

All 3 ≤ n ≤ N with σ(n) ≥ e^γ n log log n.  Unconditionally these are 26 numbers, the
largest being 5040; RH is equivalent to there being no others, ever.
"""
function robin_violations(N::Integer)
    σ = sigma_table(N)
    eγ = exp(γ_EULER)
    return [n for n in 3:N if σ[n] ≥ eγ * n * log(log(n))]
end

"""
    lagarias_holds(n) -> Bool

Lagarias (2002): RH ⇔ σ(n) ≤ Hₙ + e^{Hₙ} log Hₙ for all n ≥ 1, Hₙ the harmonic number.
"""
function lagarias_holds(n::Integer)
    H = sum(k -> 1 / k, 1:n)
    return divisor_sigma(n) ≤ H + exp(H) * log(H) + 1e-9 * n
end

"""
    lagarias_violations(N) -> Vector{Int}
"""
function lagarias_violations(N::Integer)
    σ = sigma_table(N)
    H = 0.0
    bad = Int[]
    for n in 1:N
        H += 1 / n
        σ[n] > H + exp(H) * log(H) + 1e-9 * n && push!(bad, n)
    end
    return bad
end

"""
    LI_LAMBDA1

λ₁ = 1 + γ/2 − ½ log 4π ≈ 0.0230957, the first Keiper–Li coefficient (exact).
"""
const LI_LAMBDA1 = 1 + γ_EULER / 2 - log(4π) / 2

"""
    li_coefficient(n, γs)

Keiper–Li coefficient λₙ = Σ_ρ [1 − (1 − 1/ρ)ⁿ], truncated to the zeros 1/2 ± iγ, γ ∈ `γs`.
Li (1997): RH ⇔ λₙ > 0 for all n ≥ 1. (Truncation underestimates; compare `LI_LAMBDA1`.)
"""
function li_coefficient(n::Integer, γs::AbstractVector)
    s = 0.0
    for γ in γs
        ρ = complex(0.5, γ)
        s += 2 * real(1 - (1 - 1 / ρ)^n)
    end
    return s
end

"""
    mertens_ratio_max(N) -> (x, ratio)

max_{2≤x≤N} |M(x)|/√x.  Mertens conjectured < 1 for all x > 1; Odlyzko & te Riele (1985) disproved it,
but RH ⇔ M(x) = O(x^{1/2+ε}).
"""
function mertens_ratio_max(N::Integer)
    M = mertens_table(N)
    r = abs.(M[2:N]) ./ sqrt.(2:N)
    i = argmax(r)
    return (x = i + 1, ratio = r[i])
end

"""
    zero_density_exponents(σ) -> NamedTuple

Exponents E(σ) in zero-density bounds N(σ, T) ≪ T^{E(σ)+ε}, where N(σ,T) counts zeros with
Re ρ ≥ σ, 0 < Im ρ ≤ T.  RH says N(σ,T) = 0 for σ > 1/2.

* `trivial`        : 1 (all zeros, N(T) ~ (T/2π) log T)
* `ingham`         : 3(1−σ)/(2−σ)                        (1940)
* `huxley`         : 3(1−σ)/(3σ−1)                       (1972)
* `guth_maynard`   : 30(1−σ)/13                           (2024)
* `density_hyp`    : 2(1−σ)                               (Density Hypothesis, conjectural)
"""
function zero_density_exponents(σ::Real)
    return (trivial = 1.0,
            ingham = 3(1 - σ) / (2 - σ),
            huxley = 3(1 - σ) / (3σ - 1),
            guth_maynard = 30(1 - σ) / 13,
            density_hyp = 2(1 - σ))
end
