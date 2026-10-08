# Getting started

## Installation

Riemannian is not registered. Install it from GitHub:

```julia
using Pkg
Pkg.add(url = "https://github.com/geekymode/Riemannian.jl")
Pkg.add("CairoMakie")          # optional, for plots
```

Or, from a local clone:

```julia
] dev /path/to/Riemannian.jl
```

Julia ≥ 1.10 is required, because plotting is implemented as a package extension.

## A first session

```@example start
using Riemannian
explain()
```

```@example start
explain(:trivial_zeros)
```

Evaluate ζ anywhere in the complex plane:

```@example start
zeta(2), π^2 / 6
```

```@example start
zeta(-1)          # the "−1/12"
```

```@example start
zeta(0.5 + 14.134725141734693im)   # first nontrivial zero
```

Compute zeros and verify that they lie on the critical line:

```@example start
nontrivial_zeros(5)
```

```@example start
check_zeros(500)
```

## Plotting

```@example start
using CairoMakie
set_theme!(riemann_theme())
plot_critical_line(tmax = 50)
```

Every `plot_*` function returns a `Makie.Figure`, so you can save it or modify it:

```julia
fig = plot_zeta_spiral()
save("spiral.png", fig; px_per_unit = 2)
```

To render every figure at once:

```julia
gallery("figures")             # 21 PNGs
gallery("figures"; quick = true)  # lower resolution, faster
```

## High precision

All the core ζ machinery is generic in the floating-point type:

```@example start
setprecision(BigFloat, 200) do
    γ₁ = parse(BigFloat, "14.134725141734693790457251983562470270784257115699")
    abs(zeta(complex(big"0.5", γ₁)))
end
```

## Performance notes

* `zeta(s)` costs ``O(|s|)`` operations (Euler–Maclaurin). Near ``t \sim 10^3`` it takes about 5 µs.
* `nontrivial_zeros(n)` caches its results. The first 5000 zeros take about 3 s.
* [`riemann_siegel_Z`](@ref) costs ``O(\sqrt{t})`` and is the way to go at large heights.
* The heaviest plots are `plot_heat_flow` (BigFloat integrals, about 12 s) and `plot_spacing_distribution`.
