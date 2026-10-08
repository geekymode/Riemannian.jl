# A guided tour of Riemannian — run chunk by chunk in the REPL / VS Code.
#
#   julia --project=examples -t auto
#   julia> include("examples/tour.jl")
#
using Riemannian
using CairoMakie            # activates the plotting extension
set_theme!(riemann_theme())

explain()                   # the reading list

# ════════════════════════════════════════════════════════════════════════════
# 1. Primes
# ════════════════════════════════════════════════════════════════════════════
explain(:primes)
primepi(10^6), li(1e6), riemann_R(1e6)       # 78498 vs 78627.5 vs 78527.4
chebyshev_psi(1000) / 1000                    # → 1 (PNT)
mertens(10^5), mertens_ratio_max(10^6)
plot_prime_counting()
plot_prime_gaps()

# ════════════════════════════════════════════════════════════════════════════
# 2. Euler product and special values
# ════════════════════════════════════════════════════════════════════════════
explain(:euler_product)
euler_product_partial(2, 10^4), dirichlet_partial(2, 10^4), π^2 / 6
zeta_even_exact(1), zeta_even_exact(2), zeta_even_exact(6)   # ζ(2n)/π^{2n}
plot_euler_product()

# ════════════════════════════════════════════════════════════════════════════
# 3. Analytic continuation and the functional equation
# ════════════════════════════════════════════════════════════════════════════
explain(:analytic_continuation)
s = 0.3 + 5im
zeta_borwein(s), zeta_em(s), zeta(s)          # three routes, one function
zeta(-1), zeta_negint_exact(1)                # −1/12 — of the *continued* function
explain(:functional_equation)
abs(chi_factor(0.5 + 30im))                   # = 1 on the critical line
completed_zeta(s) ≈ completed_zeta(1 - s)
plot_analytic_continuation()
plot_domain_coloring()

# ════════════════════════════════════════════════════════════════════════════
# 4. Trivial zeros
# ════════════════════════════════════════════════════════════════════════════
explain(:trivial_zeros)
trivial_zeros(5), zeta.(-2.0:-2.0:-10.0)
[zeta_negint_exact(n) for n in 0:7]
plot_zeta_real()

# ════════════════════════════════════════════════════════════════════════════
# 5. Nontrivial zeros
# ════════════════════════════════════════════════════════════════════════════
explain(:nontrivial_zeros)
γ = nontrivial_zeros(10)
hardy_Z.(γ)                                   # ≈ 0
zero_count(1000), check_zeros(1000)           # every zero up to T = 1000 is on the line
gram_point(0), gram_law_violations(300)       # Gram's law first fails at n = 126
riemann_siegel_Z(10^4), hardy_Z(10^4)         # O(√t) vs O(t) work
plot_critical_line()
plot_zeta_spiral()
plot_zero_counting()
plot_xi()
plot_modulus_surface()

# ════════════════════════════════════════════════════════════════════════════
# 6. Explicit formula — primes from zeros
# ════════════════════════════════════════════════════════════════════════════
explain(:explicit_formula)
γ = nontrivial_zeros(500)
psi_explicit(100.5, γ), chebyshev_psi(100.5)
primepi_explicit(100.5, γ), primepi(100)
plot_explicit_formula()
plot_prime_staircase()

# ════════════════════════════════════════════════════════════════════════════
# 7. Random matrices
# ════════════════════════════════════════════════════════════════════════════
explain(:random_matrices)
γ = nontrivial_zeros(5000)
close_pairs(γ; k = 5)                         # Lehmer-type near collisions
plot_spacing_distribution()
plot_pair_correlation()

# ════════════════════════════════════════════════════════════════════════════
# 8. Equivalent criteria, bounds, de Bruijn–Newman
# ════════════════════════════════════════════════════════════════════════════
explain(:equivalents)
robin_violations(10^5)                        # 26 numbers, all ≤ 5040
lagarias_violations(10^4)                     # none
li_coefficient(1, γ), LI_LAMBDA1              # truncated sum creeps up to λ₁
zero_density_exponents(0.75)
plot_robin()
plot_zero_density()
plot_heat_flow()

# ════════════════════════════════════════════════════════════════════════════
# 9. History
# ════════════════════════════════════════════════════════════════════════════
explain(:breakthroughs)
timeline()
timeline(category = :zeros)
plot_timeline()

# Everything at once, as PNGs:
# gallery("figures")
