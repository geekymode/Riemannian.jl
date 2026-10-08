```@meta
CurrentModule = Riemannian
```

# Riemannian.jl

*A computational companion to the Riemann Hypothesis.*

Riemannian is a Julia package for **studying** the Riemann Hypothesis (RH) by computing and drawing
things. It starts with counting primes, follows Euler and Riemann into the complex plane, and arrives at the
zeros of ``\zeta(s)``: the trivial ones, the nontrivial ones, how to find and count them, what they
say about primes, how they resemble eigenvalues of random matrices, and what has been proved so far.

```@example index
using Riemannian, CairoMakie
set_theme!(riemann_theme())  # hide
plot_domain_coloring(n = 360)
```

## The hypothesis in one paragraph

The Riemann zeta function is defined for ``\operatorname{Re} s > 1`` by

```math
\zeta(s) \;=\; \sum_{n=1}^{\infty} \frac{1}{n^{s}} \;=\; \prod_{p\ \text{prime}} \frac{1}{1-p^{-s}},
```

and extends to a meromorphic function on all of ``\mathbb{C}`` with a single pole at ``s = 1``.
It vanishes at ``s = -2, -4, -6, \dots`` (the *trivial zeros*). Every other zero lies in the
*critical strip* ``0 < \operatorname{Re} s < 1``.

!!! info "The Riemann Hypothesis (1859)"
    Every nontrivial zero of ``\zeta(s)`` has real part exactly ``\tfrac12``.

RH is equivalent to the statement that the primes are distributed as regularly as possible:
``\pi(x) = \operatorname{li}(x) + O(\sqrt{x}\log x)``.

## What's inside

| Chapter | What you compute | Key functions |
|:--|:--|:--|
| [1. Primes](@ref ch-primes) | ``\pi(x)``, ``\operatorname{li}(x)``, ``R(x)``, ``\psi(x)``, ``M(x)`` | [`primepi`](@ref), [`riemann_R`](@ref), [`chebyshev_psi`](@ref) |
| [2. Euler product](@ref ch-euler) | ``\zeta(2n)`` exactly, sum vs product | [`zeta_even_exact`](@ref), [`euler_product_partial`](@ref) |
| [3. Analytic continuation](@ref ch-continuation) | ``\zeta(s)`` on all of ``\mathbb{C}`` | [`zeta`](@ref), [`zeta_em`](@ref), [`zeta_borwein`](@ref), [`chi_factor`](@ref) |
| [4. Trivial zeros](@ref ch-trivial) | ``\zeta(-2n)=0``, ``\zeta(-1) = -\tfrac1{12}`` | [`trivial_zeros`](@ref), [`zeta_negint_exact`](@ref) |
| [5. Nontrivial zeros](@ref ch-nontrivial) | ``Z(t)``, Gram points, ``N(T)`` | [`hardy_Z`](@ref), [`nontrivial_zeros`](@ref), [`zero_count`](@ref) |
| [6. Explicit formulas](@ref ch-explicit) | primes rebuilt from zeros | [`psi_explicit`](@ref), [`primepi_explicit`](@ref) |
| [7. Random matrices](@ref ch-rmt) | spacings, pair correlation | [`normalized_spacings`](@ref), [`pair_correlation`](@ref) |
| [8. Equivalents](@ref ch-equivalents) | Robin, Lagarias, Li, zero density | [`robin_violations`](@ref), [`li_coefficient`](@ref) |
| [9. de Bruijn–Newman](@ref ch-dbn) | heat flow of zeros, ``\Lambda`` | [`debruijn_H`](@ref), [`heat_flow_zeros`](@ref) |
| [10. History](@ref ch-history) | 34 milestones, 1737–2024 | [`breakthroughs`](@ref), [`timeline`](@ref), [`explain`](@ref) |

Every chapter is runnable: the code blocks on these pages are executed when the documentation is
built, and every figure you see was produced by them.

## Design

* **Numerics are self-contained and generic.** `zeta`, `loggamma_c`, `xi`, … are written for any
  `AbstractFloat`, so `BigFloat` works out of the box. The only numerical dependency is
  SpecialFunctions.jl (for `expint`).
* **Plotting is optional.** CairoMakie is a *weak dependency*: `using CairoMakie` loads the
  `RiemannianCairoMakieExt` extension and all `plot_*` functions. Without it the package loads in a
  fraction of a second.
* **Learning in the REPL.** `explain()` prints a reading list, and `explain(:topic)` prints a short essay
  that ends with functions to try.
