# [1 · Primes](@id ch-primes)

```@setup ch1
using Riemannian, CairoMakie
set_theme!(riemann_theme())
```

## Counting primes

The prime-counting function ``\pi(x) = \#\{p \le x\}`` is the object of the whole story.

```@example ch1
sieve(50)
```

```@example ch1
[primepi(10^k) for k in 1:7]
```

Near ``x``, roughly one integer in ``\log x`` is prime. Gauss (around 1792) turned this observation into the
guess

```math
\pi(x) \approx \operatorname{li}(x) = \int_0^x \frac{dt}{\log t},
```

and the **Prime Number Theorem** (Hadamard and de la Vallée Poussin, 1896) states that
``\pi(x) \sim x / \log x``. Riemann proposed the better approximation

```math
R(x) = \sum_{n\ge1} \frac{\mu(n)}{n}\, \operatorname{li}\!\left(x^{1/n}\right)
     = 1 + \sum_{k\ge1} \frac{(\log x)^k}{k \cdot k!\, \zeta(k+1)} ,
```

which [`riemann_R`](@ref) evaluates with the rapidly convergent Gram series on the right.

```@example ch1
x = 10^6
(π = primepi(x), x_over_log = x / log(x), li = li(x), R = riemann_R(x))
```

```@example ch1
plot_prime_counting()
```

The right panel is the real subject of RH: the **error** ``\operatorname{li}(x) - \pi(x)``.
It is tiny compared with ``x/\log x - \pi(x)``. RH is equivalent to (von Koch, 1901)

```math
\pi(x) = \operatorname{li}(x) + O\!\left(\sqrt{x}\,\log x\right).
```

## Prime gaps

```@example ch1
ps, gs = prime_gaps(10^6)
maximum(gs), ps[argmax(gs)]
```

The average gap near ``p`` is ``\log p``. Cramér's random model predicts that the largest gaps are
about ``\log^2 p``.

```@example ch1
plot_prime_gaps()
```

## Arithmetic functions that encode primes

Several arithmetic functions have Dirichlet series that are simple expressions in ``\zeta``:

| function | definition | Dirichlet series |
|:--|:--|:--|
| Möbius ``\mu(n)`` | ``(-1)^k`` if ``n`` is a product of ``k`` distinct primes, else 0 | ``1/\zeta(s)`` |
| Liouville ``\lambda(n)`` | ``(-1)^{\Omega(n)}`` | ``\zeta(2s)/\zeta(s)`` |
| von Mangoldt ``\Lambda(n)`` | ``\log p`` if ``n = p^k``, else 0 | ``-\zeta'(s)/\zeta(s)`` |
| Euler ``\varphi(n)`` | ``\#\{k \le n : \gcd(k,n) = 1\}`` | ``\zeta(s-1)/\zeta(s)`` |
| ``\sigma_k(n)`` | ``\sum_{d \mid n} d^k`` | ``\zeta(s)\zeta(s-k)`` |

```@example ch1
[mobius(n) for n in 1:12]
```

```@example ch1
vonmangoldt(8), vonmangoldt(12), totient(36), divisor_sigma(12)
```

### Chebyshev's functions

The weighted counts

```math
\vartheta(x) = \sum_{p \le x} \log p, \qquad
\psi(x) = \sum_{n\le x} \Lambda(n) = \sum_{p^k \le x} \log p
```

are more natural than ``\pi(x)``. The PNT is equivalent to ``\psi(x) \sim x``, and
``\psi`` is what the explicit formula of [Chapter 6](@ref ch-explicit) describes exactly.

```@example ch1
chebyshev_psi(10), log(lcm(1:10))     # ψ(n) = log lcm(1,…,n)
```

```@example ch1
plot_chebyshev()
```

### The Mertens function

```math
M(x) = \sum_{n \le x} \mu(n).
```

It behaves like a random walk. RH is equivalent to ``M(x) = O(x^{1/2+\varepsilon})`` for every
``\varepsilon > 0``. The stronger *Mertens conjecture* ``|M(x)| < \sqrt{x}`` is false
(Odlyzko–te Riele 1985), although no counterexample is known explicitly.

```@example ch1
mertens(10^5), mertens_ratio_max(10^6)
```

```@example ch1
plot_mertens()
```

## Functions in this chapter

[`sieve`](@ref), [`isprime`](@ref), [`primepi`](@ref), [`primepi_table`](@ref), [`nthprime`](@ref),
[`factorize`](@ref), [`prime_gaps`](@ref), [`twin_primes`](@ref), [`li`](@ref), [`Li`](@ref),
[`riemann_R`](@ref), [`mobius`](@ref), [`mobius_table`](@ref), [`liouville`](@ref), [`mertens`](@ref),
[`mertens_table`](@ref), [`vonmangoldt`](@ref), [`chebyshev_theta`](@ref), [`chebyshev_psi`](@ref),
[`totient`](@ref), [`divisor_sigma`](@ref), [`sigma_table`](@ref).
