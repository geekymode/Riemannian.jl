# ─────────────────────────────────────────────────────────────────────────────
# Chapter 5 — Statistics of zeros & random matrices (Montgomery–Dyson, Odlyzko)
# ─────────────────────────────────────────────────────────────────────────────

"""
    unfold(γs)

Map zero heights to "unfolded" positions θ(γ)/π + 1 ≈ N(γ), which have mean spacing 1.
"""
unfold(γs::AbstractVector) = [riemann_von_mangoldt(γ) for γ in γs]

"""
    normalized_spacings(γs)

Consecutive zero gaps rescaled to mean 1.  Their distribution matches the GUE
(Gaussian Unitary Ensemble) of random Hermitian matrices — Odlyzko's famous computation.
"""
normalized_spacings(γs::AbstractVector) = diff(unfold(γs))

"""
    wigner_gue(s)

GUE Wigner surmise p(s) = (32/π²) s² e^{−4s²/π}: level repulsion ∝ s² at small gaps.
"""
wigner_gue(s::Real) = 32 / π^2 * s^2 * exp(-4s^2 / π)

"""
    wigner_goe(s)

GOE (real symmetric) Wigner surmise p(s) = (π/2) s e^{−πs²/4}.
"""
wigner_goe(s::Real) = π / 2 * s * exp(-π * s^2 / 4)

"""
    poisson_spacing(s)

Spacing density e^{−s} of uncorrelated (Poisson) points — what independent random zeros would give.
"""
poisson_spacing(s::Real) = exp(-s)

"""
    montgomery_pair_correlation(u)

Montgomery's (1973) pair-correlation density 1 − (sin πu / πu)², identified by Dyson
as the GUE eigenvalue pair correlation.
"""
montgomery_pair_correlation(u::Real) = u == 0 ? 0.0 : 1 - (sinpi(u) / (π * u))^2

"""
    pair_correlation(γs; umax = 3, nbins = 60) -> (centers, density)

Empirical pair correlation of unfolded zeros: histogram of all differences uⱼ − uᵢ ∈ (0, umax],
normalised so that uncorrelated points would give density 1.
"""
function pair_correlation(γs::AbstractVector; umax::Real = 3, nbins::Integer = 60)
    u = unfold(sort(γs))
    edges = range(0, umax; length = nbins + 1)
    counts = zeros(Int, nbins)
    n = length(u)
    for i in 1:n
        for j in i+1:n          # separate loop: `break` must only leave the inner one
            d = u[j] - u[i]
            d > umax && break
            b = min(nbins, floor(Int, d / step(edges)) + 1)
            counts[b] += 1
        end
    end
    centers = (edges[1:end-1] .+ edges[2:end]) ./ 2
    return collect(centers), counts ./ (n * step(edges))
end

"""
    gue_eigenvalues(N; rng)

Eigenvalues of an N×N GUE matrix H = (A + A†)/2, A with iid standard complex Gaussian entries.
"""
function gue_eigenvalues(N::Integer; rng::AbstractRNG = Random.default_rng())
    A = randn(rng, ComplexF64, N, N)
    return eigvals(Hermitian((A + A') / 2))
end

"""
    gue_spacings(N = 400; samples = 20, rng) -> Vector{Float64}

Unfolded nearest-neighbour spacings from the bulk (central half) of GUE spectra,
for comparison with [`normalized_spacings`](@ref) of zeta zeros.
"""
function gue_spacings(N::Integer = 400; samples::Integer = 20,
                      rng::AbstractRNG = Random.default_rng())
    R = sqrt(2N)                                     # semicircle radius for this scaling
    count(x) = N * (0.5 + (x * sqrt(R^2 - x^2) / R^2 + asin(x / R)) / π)
    out = Float64[]
    for _ in 1:samples
        λ = gue_eigenvalues(N; rng)
        bulk = λ[N÷4:3N÷4]
        append!(out, diff(count.(bulk)))
    end
    return out
end
