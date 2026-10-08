# [2 · Euler product & ζ at integers](@id ch-euler)

```@setup ch2
using Riemannian, CairoMakie
set_theme!(riemann_theme())
```

## The Euler product

For ``\operatorname{Re} s > 1``,

```math
\zeta(s) = \sum_{n=1}^\infty n^{-s} = \prod_{p} \left(1 - p^{-s}\right)^{-1}.
```

*Why it holds:* expand each factor as a geometric series ``1 + p^{-s} + p^{-2s} + \cdots`` and multiply
the series together. By unique factorisation, every ``n^{-s}`` appears exactly once. The Euler product
is the fundamental theorem of arithmetic written as an identity between functions.

It has three immediate consequences:

1. ``\zeta(s) \ne 0`` for ``\operatorname{Re} s > 1``, because a convergent product of nonzero factors is nonzero.
2. ``\zeta(s) \to \infty`` as ``s \to 1^+``, so there are infinitely many primes. In fact
   ``\sum_p 1/p = \infty`` (Euler, 1737).
3. ``\log \zeta(s) = \sum_p \sum_k \frac{p^{-ks}}{k}``: zeta *is* a generating function for primes.

```@example ch2
euler_product_partial(2, 1000), dirichlet_partial(2, 1000), zeta(2)
```

The two kinds of truncation converge at different rates. Summing over ``n \le N`` leaves an error of
about ``1/N``, while multiplying over ``p \le N`` leaves an error that drops only when ``N`` passes a prime:

```@example ch2
plot_euler_product()
```

## ζ at even integers: the Basel problem

Euler showed in 1735 that ``\zeta(2) = \pi^2/6``, and more generally that

```math
\zeta(2n) = \frac{(-1)^{n+1} B_{2n} (2\pi)^{2n}}{2\,(2n)!},
```

where ``B_k`` are the Bernoulli numbers (``\frac{t}{e^t-1} = \sum B_k \frac{t^k}{k!}``).
[`zeta_even_exact`](@ref) returns the rational number ``\zeta(2n)/\pi^{2n}`` exactly:

```@example ch2
[n => zeta_even_exact(n) for n in 1:6]
```

```@example ch2
Float64(zeta_even_exact(3)) * π^6, zeta(6)
```

```@example ch2
[bernoulli(k) for k in 0:12]
```

!!! note "Odd integers are mysterious"
    No closed form is known for ``\zeta(3), \zeta(5), \dots``. Apéry (1978) proved that
    ``\zeta(3)`` is irrational. Ball–Rivoal (2001) showed that infinitely many ``\zeta(2n+1)`` are irrational.

```@example ch2
zeta(3)    # Apéry's constant
```

## Related: the Dirichlet eta function

The alternating version

```math
\eta(s) = \sum_{n\ge1} \frac{(-1)^{n-1}}{n^s} = \left(1 - 2^{1-s}\right)\zeta(s)
```

converges for ``\operatorname{Re} s > 0``. This makes it the first tool for continuing ``\zeta`` beyond
``\operatorname{Re} s > 1``; see the [next chapter](@ref ch-continuation).

```@example ch2
dirichlet_eta(1), log(2)
```

## Functions in this chapter

[`bernoulli`](@ref), [`zeta_even_exact`](@ref), [`dirichlet_partial`](@ref),
[`euler_product_partial`](@ref), [`eta_partial`](@ref), [`dirichlet_eta`](@ref).
