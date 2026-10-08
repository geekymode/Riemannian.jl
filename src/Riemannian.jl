"""
    Riemannian

A study kit for the Riemann Hypothesis: from primes and the Euler product, through
analytic continuation and the functional equation, to trivial and nontrivial zeros,
explicit formulas, random-matrix statistics, equivalent criteria and recent results.

The numerical core has no plotting dependency. Loading `CairoMakie` activates the
plotting extension (`plot_*` functions and [`gallery`](@ref)).

Start with `explain()` for the list of topics, or `explain(:trivial_zeros)` etc.
"""
module Riemannian

using LinearAlgebra
using Printf
using Random
using SpecialFunctions: expint, expinti

const γ_EULER = Base.MathConstants.eulergamma

include("primes.jl")
include("arithmetic.jl")
include("special.jl")
include("zeta.jl")
include("zeros.jl")
include("explicit.jl")
include("statistics.jl")
include("criteria.jl")
include("debruijn_newman.jl")
include("history.jl")
include("plotting_api.jl")

# 1. Primes and arithmetic functions
export sieve, isprime, primepi, primepi_table, nthprime, factorize, prime_gaps, twin_primes,
       mobius, mobius_table, liouville, mertens, mertens_table, vonmangoldt,
       chebyshev_theta, chebyshev_psi, totient, divisor_sigma, sigma_table,
       li, Li, riemann_R

# 2. Zeta function: series, products, continuation
export bernoulli, zeta_even_exact, zeta_negint_exact,
       dirichlet_partial, euler_product_partial, eta_partial,
       dirichlet_eta, zeta_borwein, zeta_em, zeta, chi_factor, completed_zeta, xi, Xi, partial_sums,
       loggamma_c

# 3. Zeros
export trivial_zeros, KNOWN_ZEROS, riemann_siegel_theta, hardy_Z, riemann_siegel_Z,
       gram_point, gram_points, gram_law_violations, argzeta_S, zero_count,
       riemann_von_mangoldt, zeros_between, nontrivial_zeros, check_zeros, close_pairs,
       RH_VERIFIED_HEIGHT, zero_free_boundary

# 4. Explicit formulas
export psi_explicit, riemann_J, riemann_J_explicit, primepi_explicit

# 5. Statistics of zeros
export unfold, normalized_spacings, wigner_gue, wigner_goe, poisson_spacing,
       montgomery_pair_correlation, pair_correlation, gue_eigenvalues, gue_spacings

# 6. Equivalent criteria / bounds
export robin_ratio, robin_violations, lagarias_holds, lagarias_violations,
       li_coefficient, LI_LAMBDA1, mertens_ratio_max, zero_density_exponents

# 7. de Bruijn–Newman
export debruijn_Phi, debruijn_H, heat_flow_zeros

# 8. History & narrative
export breakthroughs, timeline, explain

# 9. Plotting (implemented in the CairoMakie extension)
export riemann_theme, plot_prime_counting, plot_prime_gaps, plot_chebyshev, plot_mertens,
       plot_euler_product, plot_zeta_real, plot_analytic_continuation, plot_domain_coloring,
       plot_modulus_surface, plot_critical_line, plot_zeta_spiral, plot_zero_counting,
       plot_explicit_formula, plot_prime_staircase, plot_spacing_distribution,
       plot_pair_correlation, plot_xi, plot_robin, plot_zero_density, plot_heat_flow,
       plot_timeline, plot_critical_strip, plot_strip_schematic, plot_strip_width,
       plot_partial_sum_spiral, plot_zeta_near_origin, record_zeta_spiral, gallery

end # module
