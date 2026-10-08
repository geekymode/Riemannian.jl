# [8 · Equivalent criteria & bounds](@id ch-equivalents)

```@setup ch8
using Riemannian, CairoMakie
set_theme!(riemann_theme())
```

RH has hundreds of equivalent formulations. Several of them are elementary statements about integers,
and Riemannian can test those numerically. Of course, no finite computation can prove them.

## von Koch (1901)

```math
\text{RH} \iff \pi(x) = \operatorname{li}(x) + O\!\left(\sqrt x \log x\right).
```

See the error panel in [Chapter 1](@ref ch-primes).

## Mertens function

```math
\text{RH} \iff M(x) = O\!\left(x^{1/2+\varepsilon}\right) \text{ for every } \varepsilon > 0.
```

```@example ch8
mertens_ratio_max(10^6)
```

## Robin's criterion (1984)

```math
\text{RH} \iff \sigma(n) < e^{\gamma}\, n \log\log n \quad \text{for all } n > 5040.
```

Unconditionally, the inequality fails for exactly 26 values ``3 \le n \le 5040`` (27 if you also count ``n = 2``,
where ``\log\log 2 < 0``):

```@example ch8
robin_violations(10^5)
```

```@example ch8
plot_robin()
```

Akbary and Friggstad showed that if Robin's inequality fails at all, it fails at a *superabundant* number.
Computations by Briggs and by Morrill–Platt rule out every candidate up to an astronomically large bound.

## Lagarias' criterion (2002)

With ``H_n = 1 + \frac12 + \dots + \frac1n``,

```math
\text{RH} \iff \sigma(n) \le H_n + e^{H_n}\log H_n \quad \text{for all } n \ge 1.
```

This version needs no exceptions and no limits:

```@example ch8
lagarias_violations(10^4)
```

## Li's criterion (1997)

Define the Keiper–Li coefficients

```math
\lambda_n = \sum_\rho \left[1 - \left(1 - \frac1\rho\right)^n\right].
```

Then RH ``\iff \lambda_n > 0`` for all ``n \ge 1``. The first one is known in closed form,
``\lambda_1 = 1 + \frac\gamma2 - \frac12\log 4\pi \approx 0.0230957``. The truncated sum creeps up
towards it:

```@example ch8
[(k, li_coefficient(1, nontrivial_zeros(k))) for k in (10, 100, 1000, 5000)]
```

```@example ch8
LI_LAMBDA1
```

## Zero-density estimates

RH says that ``N(\sigma, T) = \#\{\rho : \beta \ge \sigma, 0 < \gamma \le T\}`` vanishes for every ``\sigma > \frac12``.
Unconditionally we have bounds of the form

```math
N(\sigma, T) \ll T^{A(\sigma)(1-\sigma) + \varepsilon}.
```

| result | exponent ``A(\sigma)(1-\sigma)`` |
|:--|:--|
| Ingham (1940) | ``\frac{3(1-\sigma)}{2-\sigma}`` |
| Huxley (1972) | ``\frac{3(1-\sigma)}{3\sigma - 1}`` |
| **Guth–Maynard (2024)** | ``\frac{30}{13}(1-\sigma)`` |
| Density Hypothesis (conjecture) | ``2(1-\sigma)`` |

The Guth–Maynard bound was the first substantial improvement in the critical range near ``\sigma = \frac34``
since Ingham. Its consequences include the prime number theorem in short intervals
``[x, x + x^{17/30 + \varepsilon}]``, improving Huxley's exponent ``7/12``.

```@example ch8
zero_density_exponents(0.75)
```

```@example ch8
plot_zero_density()
```

## Proportion of zeros on the line

Write ``\kappa`` for the proportion of nontrivial zeros on the line that have been proved.

| | ``\kappa >`` |
|:--|:--|
| Hardy (1914) | infinitely many |
| Selberg (1942) | a positive proportion |
| Levinson (1974) | ``1/3`` |
| Conrey (1989) | ``2/5`` |
| Bui–Conrey–Young (2011) | ``0.4105`` |
| Pratt–Robles–Zaharescu–Zeindler (2020) | ``5/12 \approx 0.4167`` |

## Functions in this chapter

[`robin_ratio`](@ref), [`robin_violations`](@ref), [`lagarias_holds`](@ref),
[`lagarias_violations`](@ref), [`li_coefficient`](@ref), [`LI_LAMBDA1`](@ref),
[`mertens_ratio_max`](@ref), [`zero_density_exponents`](@ref).
