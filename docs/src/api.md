# API reference

```@meta
CurrentModule = Riemannian
```

```@docs
Riemannian
```

## Primes

```@docs
sieve
isprime
primepi
primepi_table
nthprime
factorize
prime_gaps
twin_primes
li
Li
riemann_R
```

## Arithmetic functions

```@docs
mobius
mobius_table
liouville
mertens
mertens_table
vonmangoldt
chebyshev_theta
chebyshev_psi
totient
divisor_sigma
sigma_table
```

## Exact values, series and products

```@docs
bernoulli
zeta_even_exact
zeta_negint_exact
dirichlet_partial
eta_partial
euler_product_partial
```

## The zeta function and relatives

```@docs
zeta
zeta_em
zeta_borwein
dirichlet_eta
chi_factor
completed_zeta
xi
Xi
loggamma_c
partial_sums
```

## Zeros

```@docs
trivial_zeros
KNOWN_ZEROS
riemann_siegel_theta
hardy_Z
riemann_siegel_Z
gram_point
gram_points
gram_law_violations
argzeta_S
riemann_von_mangoldt
zero_count
zeros_between
nontrivial_zeros
check_zeros
close_pairs
RH_VERIFIED_HEIGHT
zero_free_boundary
```

## Explicit formulas

```@docs
psi_explicit
riemann_J
riemann_J_explicit
primepi_explicit
```

## Statistics of zeros

```@docs
unfold
normalized_spacings
wigner_gue
wigner_goe
poisson_spacing
montgomery_pair_correlation
pair_correlation
gue_eigenvalues
gue_spacings
```

## Criteria and bounds

```@docs
robin_ratio
robin_violations
lagarias_holds
lagarias_violations
LI_LAMBDA1
li_coefficient
mertens_ratio_max
zero_density_exponents
```

## de Bruijn–Newman

```@docs
debruijn_Phi
debruijn_H
heat_flow_zeros
```

## History and narrative

```@docs
breakthroughs
timeline
explain
```

## [Plotting](@id api-plotting)

These functions are defined once `CairoMakie` is loaded.

```@docs
riemann_theme
plot_prime_counting
plot_prime_gaps
plot_chebyshev
plot_mertens
plot_euler_product
plot_zeta_real
plot_analytic_continuation
plot_domain_coloring
plot_modulus_surface
plot_critical_line
plot_zeta_spiral
plot_zero_counting
plot_xi
plot_explicit_formula
plot_prime_staircase
plot_spacing_distribution
plot_pair_correlation
plot_robin
plot_zero_density
plot_heat_flow
plot_timeline
plot_critical_strip
plot_strip_schematic
plot_strip_width
plot_partial_sum_spiral
plot_zeta_near_origin
record_zeta_spiral
gallery
```

## Index

```@index
```
