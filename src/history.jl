# ─────────────────────────────────────────────────────────────────────────────
# Chapter 8 — History, milestones and narrative explanations
# ─────────────────────────────────────────────────────────────────────────────

"""
    breakthroughs() -> Vector{NamedTuple}

Milestones on the road to (and around) the Riemann Hypothesis, each with
`year`, `who`, `what`, `category` (:foundations, :zeros, :computation, :equivalents,
:bounds, :analogy).
"""
breakthroughs() = [
    (year = 1737, who = "Euler", category = :foundations,
     what = "Euler product Σ n⁻ˢ = Π (1 − p⁻ˢ)⁻¹; divergence of Σ 1/p"),
    (year = 1792, who = "Gauss", category = :foundations,
     what = "Conjectures π(x) ≈ li(x) (density of primes ≈ 1/log x)"),
    (year = 1850, who = "Chebyshev", category = :bounds,
     what = "π(x) is within ~10% of x/log x; introduces ϑ(x) and ψ(x)"),
    (year = 1859, who = "Riemann", category = :foundations,
     what = "ζ(s) on ℂ: continuation, functional equation, explicit formula, the Hypothesis"),
    (year = 1895, who = "von Mangoldt", category = :foundations,
     what = "Proves the explicit formula for ψ(x) (and, in 1905, N(T) ~ (T/2π) log(T/2πe))"),
    (year = 1896, who = "Hadamard; de la Vallée Poussin", category = :bounds,
     what = "Prime Number Theorem, via ζ(1 + it) ≠ 0"),
    (year = 1901, who = "von Koch", category = :equivalents,
     what = "RH ⇔ π(x) = li(x) + O(√x log x)"),
    (year = 1903, who = "Gram", category = :computation,
     what = "First 15 zeros computed; Gram points and 'Gram's law'"),
    (year = 1914, who = "Hardy", category = :zeros,
     what = "Infinitely many zeros lie on the critical line"),
    (year = 1914, who = "Littlewood", category = :bounds,
     what = "π(x) − li(x) changes sign infinitely often"),
    (year = 1932, who = "Siegel", category = :computation,
     what = "Riemann–Siegel formula, recovered from Riemann's unpublished notes"),
    (year = 1940, who = "Ingham", category = :bounds,
     what = "Zero-density estimate N(σ,T) ≪ T^(3(1−σ)/(2−σ)+ε)"),
    (year = 1942, who = "Selberg", category = :zeros,
     what = "A positive proportion of zeros lie on the critical line"),
    (year = 1948, who = "Weil", category = :analogy,
     what = "RH for curves over finite fields proved"),
    (year = 1953, who = "Turing", category = :computation,
     what = "Turing's method to certify zero counts; early machine computation"),
    (year = 1972, who = "Huxley", category = :bounds,
     what = "Zero-density N(σ,T) ≪ T^(12(1−σ)/5+ε); primes in short intervals x^(7/12)"),
    (year = 1973, who = "Montgomery; Dyson", category = :zeros,
     what = "Pair correlation 1 − (sin πu/πu)²: zeros look like GUE eigenvalues"),
    (year = 1974, who = "Levinson", category = :zeros,
     what = "More than 1/3 of zeros lie on the critical line"),
    (year = 1974, who = "Deligne", category = :analogy,
     what = "Weil conjectures: RH for varieties over finite fields"),
    (year = 1976, who = "Newman", category = :equivalents,
     what = "Conjectures Λ ≥ 0 for the de Bruijn–Newman constant (RH ⇔ Λ ≤ 0)"),
    (year = 1984, who = "Robin", category = :equivalents,
     what = "RH ⇔ σ(n) < e^γ n log log n for all n > 5040"),
    (year = 1985, who = "Odlyzko; te Riele", category = :equivalents,
     what = "Disprove the Mertens conjecture |M(x)| < √x"),
    (year = 1987, who = "Odlyzko", category = :computation,
     what = "Zeros near the 10²⁰-th match GUE spacing statistics to high accuracy"),
    (year = 1989, who = "Conrey", category = :zeros,
     what = "More than 2/5 of zeros lie on the critical line"),
    (year = 1997, who = "Li", category = :equivalents,
     what = "Li's criterion: RH ⇔ λₙ > 0 for all n"),
    (year = 2000, who = "Clay Mathematics Institute", category = :foundations,
     what = "RH named a Millennium Prize Problem"),
    (year = 2002, who = "Lagarias", category = :equivalents,
     what = "RH ⇔ σ(n) ≤ Hₙ + e^(Hₙ) log Hₙ for all n"),
    (year = 2004, who = "Gourdon", category = :computation,
     what = "First 10¹³ zeros verified on the critical line"),
    (year = 2011, who = "Bui; Conrey; Young", category = :zeros,
     what = "More than 41% of zeros on the critical line"),
    (year = 2019, who = "Polymath15", category = :equivalents,
     what = "Λ ≤ 0.22"),
    (year = 2020, who = "Rodgers; Tao", category = :equivalents,
     what = "Newman's conjecture Λ ≥ 0 proved: RH, if true, is 'barely' true"),
    (year = 2020, who = "Pratt; Robles; Zaharescu; Zeindler", category = :zeros,
     what = "More than 5/12 of zeros on the critical line"),
    (year = 2021, who = "Platt; Trudgian", category = :computation,
     what = "RH verified up to height 3·10¹²; Λ ≤ 0.2"),
    (year = 2024, who = "Guth; Maynard", category = :bounds,
     what = "N(σ,T) ≪ T^(30(1−σ)/13+ε): first big improvement on Ingham's 1940 bound"),
]

"""
    timeline(; category = nothing, io = stdout)

Print the milestones from [`breakthroughs`](@ref), optionally filtered by category.
"""
function timeline(; category::Union{Nothing,Symbol} = nothing, io::IO = stdout)
    for b in breakthroughs()
        (category === nothing || b.category == category) || continue
        @printf(io, "%5d  %-34s %s\n", b.year, b.who, b.what)
    end
end

const _TOPICS = Dict{Symbol,String}(
:primes => """
PRIMES AND THEIR DENSITY
  The primes thin out: near x roughly one integer in log x is prime. Gauss guessed
  π(x) ≈ li(x) = ∫ dt/log t; the Prime Number Theorem (1896) says π(x) ~ x/log x.
  The real question is the ERROR π(x) − li(x). Riemann's insight: that error is a sum
  of waves, one per zero of ζ. RH ⇔ the error is O(√x log x).
  Try: primepi(10^6), li(1e6), riemann_R(1e6), chebyshev_psi(1000), mertens(1000)
       plot_prime_counting(), plot_chebyshev(), plot_prime_gaps()
""",
:euler_product => """
THE EULER PRODUCT
  For Re s > 1:   ζ(s) = Σ_{n≥1} n⁻ˢ = Π_p (1 − p⁻ˢ)⁻¹.
  Expanding each factor as a geometric series and multiplying out gives every n⁻ˢ exactly
  once — this is unique factorisation written analytically. Consequences: ζ(s) ≠ 0 for
  Re s > 1; ζ(1) = ∞ forces infinitely many primes; log ζ(s) ≈ Σ_p p⁻ˢ.
  Try: euler_product_partial(2, 1000), dirichlet_partial(2, 1000), zeta(2) - π^2/6
       plot_euler_product()
""",
:zeta => """
THE ZETA FUNCTION AT INTEGERS
  Euler: ζ(2) = π²/6, and generally ζ(2n) = rational × π²ⁿ via Bernoulli numbers.
  ζ(3), ζ(5), … are mysterious (Apéry 1978: ζ(3) is irrational).
  Try: zeta_even_exact(1), zeta_even_exact(3), bernoulli(12), zeta(3)
""",
:analytic_continuation => """
ANALYTIC CONTINUATION
  Σ n⁻ˢ diverges for Re s ≤ 1, but the function it defines extends uniquely to a
  meromorphic function on ℂ with a single simple pole at s = 1 (residue 1).
  Three routes, all in this package:
    1. Alternating series: ζ(s) = η(s)/(1 − 2^{1−s}), η converges for Re s > 0
       (zeta_borwein, eta_partial).
    2. Euler–Maclaurin: replace the tail of the sum by an integral plus Bernoulli
       corrections; every term is analytic for s ≠ 1 (zeta_em).
    3. Functional equation: reflect Re s < 1/2 to Re s > 1/2 (zeta, chi_factor).
  "1 + 2 + 3 + ⋯ = −1/12" is shorthand for ζ(−1) = −1/12 of the CONTINUED function.
  Try: zeta(-1), zeta_negint_exact(1), dirichlet_partial(0.5, 10^4), zeta(0.5)
       plot_analytic_continuation(), plot_zeta_real()
""",
:functional_equation => """
THE FUNCTIONAL EQUATION
  ζ(s) = χ(s) ζ(1−s),  χ(s) = 2ˢ π^{s−1} sin(πs/2) Γ(1−s).
  Symmetric form: Λ(s) = π^{−s/2} Γ(s/2) ζ(s) satisfies Λ(s) = Λ(1−s).
  So ζ's behaviour is mirrored across the critical line Re s = 1/2. On that line
  |χ| = 1, so |ζ(1/2+it)| = |ζ(1/2−it)|, and Z(t) = e^{iθ(t)} ζ(1/2+it) is real.
  Try: chi_factor(0.5 + 10im) |> abs, completed_zeta(0.3+2im) - completed_zeta(0.7-2im)
""",
:trivial_zeros => """
TRIVIAL ZEROS: s = −2, −4, −6, …
  In ζ(s) = χ(s) ζ(1−s), for s = −2n (n ≥ 1) the factor sin(πs/2) = sin(−nπ) = 0, while
  Γ(1−s) and ζ(1−s) are finite and nonzero. Hence ζ(−2n) = 0. At s = 0, −1, −3, … the
  sine is nonzero, and ζ(−n) = (−1)ⁿ Bₙ₊₁/(n+1): ζ(0) = −1/2, ζ(−1) = −1/12, ζ(−3) = 1/120.
  In the completed Λ(s) these zeros are cancelled by the poles of Γ(s/2): they are
  "trivial" because they come from the Gamma factor, not from the primes.
  For s < 0, |ζ(s)| grows factorially between trivial zeros (Γ(1−s) factor).
  Try: trivial_zeros(5), zeta(-4.0), zeta(-3.0), zeta_negint_exact(3)
       plot_zeta_real()
""",
:nontrivial_zeros => """
NONTRIVIAL ZEROS AND THE CRITICAL STRIP
  All other zeros lie in the critical strip 0 < Re s < 1 (the Euler product rules out
  Re s > 1, Hadamard–de la Vallée Poussin rule out Re s = 1, the functional equation
  mirrors to Re s < 0). They are symmetric under s ↦ 1−s and s ↦ s̄.
  THE RIEMANN HYPOTHESIS: every nontrivial zero has Re s = 1/2.
  First zeros: 1/2 ± 14.1347i, 1/2 ± 21.0220i, 1/2 ± 25.0109i, …
  Counting: N(T) = θ(T)/π + 1 + S(T) ≈ (T/2π) log(T/2πe) + 7/8.
  Verification: count sign changes of the real function Z(t) and compare with N(T)
  from the argument principle. If they agree, every zero up to T is on the line.
  Try: nontrivial_zeros(10), hardy_Z(14.134725), zero_count(100), check_zeros(1000)
       gram_law_violations(200), plot_critical_line(), plot_domain_coloring(),
       plot_zeta_spiral(), plot_zero_counting(), plot_xi()
""",
:explicit_formula => """
THE EXPLICIT FORMULA: PRIMES = SMOOTH TERM − WAVES FROM ZEROS
  ψ(x) = x − Σ_ρ x^ρ/ρ − log 2π − ½ log(1 − x⁻²).
  Each pair ρ = 1/2 ± iγ contributes ≈ 2√x cos(γ log x − arg ρ)/|ρ|: a wave of frequency γ
  in log x. Summing more zeros sharpens the curve into the prime staircase.
  A zero off the line with Re ρ = β > 1/2 would produce a wave of size x^β, which is why
  RH ⇔ "the primes are as regular as possible".
  Try: psi_explicit(100.5, nontrivial_zeros(100)), primepi_explicit(100.5, nontrivial_zeros(200))
       plot_explicit_formula(), plot_prime_staircase()
""",
:random_matrices => """
ZEROS AND RANDOM MATRICES
  Montgomery (1973) computed (assuming RH) the pair correlation of zeros:
  1 − (sin πu / πu)². Dyson recognised it as the pair correlation of eigenvalues of large
  random Hermitian matrices (GUE). Odlyzko's computations near the 10²⁰-th zero confirm
  the spacing distribution with striking accuracy. Zeros repel: tiny gaps are rare (∝ s²).
  Hilbert–Pólya dream: zeros = eigenvalues of a self-adjoint operator, which would force
  them onto the line.
  Try: normalized_spacings(nontrivial_zeros(2000)), gue_spacings(), close_pairs(...)
       plot_spacing_distribution(), plot_pair_correlation()
""",
:equivalents => """
STATEMENTS EQUIVALENT TO RH
  • von Koch:  π(x) = li(x) + O(√x log x).
  • Mertens function: M(x) = O(x^{1/2+ε}) (but |M(x)| < √x is FALSE, Odlyzko–te Riele).
  • Robin:     σ(n) < e^γ n log log n for all n > 5040.
  • Lagarias:  σ(n) ≤ Hₙ + e^{Hₙ} log Hₙ for all n ≥ 1.
  • Li:        λₙ = Σ_ρ [1 − (1 − 1/ρ)ⁿ] > 0 for all n.
  • de Bruijn–Newman: Λ ≤ 0 (and Rodgers–Tao proved Λ ≥ 0, so RH ⇔ Λ = 0).
  Try: robin_violations(10^5), lagarias_violations(10^4), mertens_ratio_max(10^6),
       li_coefficient(1, nontrivial_zeros(2000)) vs LI_LAMBDA1
       plot_robin(), plot_mertens(), plot_heat_flow()
""",
:breakthroughs => """
RECENT AND LANDMARK PROGRESS
  • Proportion of zeros on the line: Selberg (>0), Levinson (>1/3), Conrey (>2/5),
    Bui–Conrey–Young (>41%), Pratt–Robles–Zaharescu–Zeindler (>5/12).
  • Computation: Platt–Trudgian (2021) — RH holds up to height 3·10¹².
  • de Bruijn–Newman: Rodgers–Tao (2020) Λ ≥ 0; Polymath15 Λ ≤ 0.22; Platt–Trudgian Λ ≤ 0.2.
  • Zero density: Guth–Maynard (2024) N(σ,T) ≪ T^{30(1−σ)/13+ε}, improving Ingham (1940);
    consequences include the PNT in short intervals [x, x + x^{17/30+ε}].
  • Function-field analogue proved (Weil, Deligne).
  Try: timeline(), timeline(category = :zeros), plot_timeline(), plot_zero_density()
""",
)

"""
    explain()
    explain(topic::Symbol)

Short guided explanations. `explain()` lists topics in reading order.
"""
function explain(io::IO = stdout)
    println(io, "Topics (in reading order):")
    for t in (:primes, :euler_product, :zeta, :analytic_continuation, :functional_equation,
              :trivial_zeros, :nontrivial_zeros, :explicit_formula, :random_matrices,
              :equivalents, :breakthroughs)
        println(io, "  explain(:", t, ")")
    end
end
function explain(topic::Symbol; io::IO = stdout)
    haskey(_TOPICS, topic) || throw(ArgumentError("unknown topic $topic; see explain()"))
    print(io, _TOPICS[topic])
end
