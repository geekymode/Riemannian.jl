# Gallery

```@setup gal
using Riemannian, CairoMakie
set_theme!(riemann_theme())
```

Every figure in Riemannian, with the call that produced it. All plotting functions take keyword
arguments; see the [API reference](@ref api-plotting). Continuous colour scales use **viridis**.

## Primes

```@example gal
plot_prime_counting()
```

```@example gal
plot_prime_gaps()
```

```@example gal
plot_chebyshev()
```

```@example gal
plot_mertens()
```

## ζ on the real line and analytic continuation

```@example gal
plot_euler_product()
```

```@example gal
plot_zeta_real()
```

```@example gal
plot_analytic_continuation()
```

## ζ on the complex plane

```@example gal
plot_domain_coloring(re = (-30, 10), im = (-5, 60), n = 500)
```

```@example gal
plot_modulus_surface(re = (-0.5, 1.5), im = (10, 50))
```

## The critical strip

```@example gal
plot_strip_schematic()
```

```@example gal
plot_critical_strip(tmax = 150)
```

```@example gal
plot_strip_width()
```

## The critical line

```@example gal
plot_critical_line(tmax = 100)
```

```@example gal
plot_zeta_spiral(tmax = 35)
```

```@example gal
plot_zero_counting(T = 200)
```

```@example gal
plot_xi(tmax = 40)
```

## Explicit formulas

```@example gal
plot_explicit_formula(xmax = 100, ks = (10, 50, 300))
```

```@example gal
plot_prime_staircase(xmax = 100, ks = (20, 100, 400), npts = 500)
```

## Statistics

```@example gal
plot_spacing_distribution(nzeros = 5000)
```

```@example gal
plot_pair_correlation(nzeros = 5000)
```

## Criteria and bounds

```@example gal
plot_robin(N = 2 * 10^5)
```

```@example gal
plot_zero_density()
```

```@example gal
plot_heat_flow(ts = range(-8, 1; length = 19), zmax = 80)
```

## History

```@example gal
plot_timeline()
```

## Saving everything

```julia
gallery("figures"; px_per_unit = 2)
```
