# [5 · Nontrivial zeros](@id ch-nontrivial)

```@setup ch5
using Riemannian, CairoMakie
set_theme!(riemann_theme())
```

## The critical strip

The nontrivial zeros ``\rho = \beta + i\gamma`` all lie in ``0 < \beta < 1``:

* The Euler product rules out ``\beta > 1``.
* Hadamard and de la Vallée Poussin (1896) ruled out ``\beta = 1``. This fact *is* the Prime Number Theorem.
* The functional equation reflects ``\beta < 0`` to ``\beta > 1``.

They are symmetric under ``\rho \mapsto 1 - \rho`` and under ``\rho \mapsto \bar\rho``. A zero off the critical line
would therefore come with three partners, ``\bar\rho``, ``1-\rho`` and ``1-\bar\rho``.

!!! info "The Riemann Hypothesis"
    ``\beta = \tfrac12`` for every nontrivial zero.

```@example ch5
plot_strip_schematic()
```

Turning the picture on its side, here is the band itself over the landscape ``\log|\zeta(\sigma + it)|``.
Every zero is a dark well, and every well sits on the line ``\sigma = \tfrac12``:

```@example ch5
plot_critical_strip()
```

```@example ch5
γ = nontrivial_zeros(10)
```

```@example ch5
γ .- collect(KNOWN_ZEROS)       # agreement with Odlyzko's tables
```

## Hardy's Z-function

How can a computation *prove* that a zero lies exactly on the line? Define the Riemann–Siegel theta
function and Hardy's ``Z``:

```math
\theta(t) = \arg\Gamma\!\left(\tfrac14 + \tfrac{it}{2}\right) - \frac{t}{2}\log\pi,
\qquad
Z(t) = e^{i\theta(t)}\,\zeta\!\left(\tfrac12 + it\right).
```

The functional equation forces ``Z(t)`` to be **real** for real ``t``, with ``|Z(t)| = |\zeta(\tfrac12+it)|``.
A real continuous function that changes sign must vanish, so each sign change of ``Z`` is a zero
*exactly* on the critical line. No precision argument is needed beyond determining the sign.

```@example ch5
hardy_Z.(γ[1:3])
```

```@example ch5
plot_critical_line()
```

Asymptotically ``\theta(t) \approx \frac t2 \log\frac{t}{2\pi e} - \frac\pi8 + \frac{1}{48t}``:

```@example ch5
t = 1000.0
riemann_siegel_theta(t), t/2 * log(t/(2π*ℯ)) - π/8 + 1/(48t)
```

### The spiral picture

As ``t`` increases, ``\zeta(\tfrac12 + it)`` traces loops that pass **through the origin** at every zero.
Off the line, at ``\sigma = 0.6``, the curve misses the origin, at least for the heights shown:

```@example ch5
plot_zeta_spiral()
```

### Spirals near zero

There is a second spiral picture, popularised by Quanta Magazine's explainer. Fix ``s = \tfrac12 + it`` and draw
the terms ``1, 2^{-s}, 3^{-s}, \dots`` head to tail, as a chain of vectors. The ``n``-th vector has length
``n^{-1/2}`` and turns by ``-t\log n``. Once ``n \gtrsim t/2\pi`` consecutive vectors point in nearly the same
direction, and the chain becomes a smooth spiral.

On the critical line the series diverges, so the raw chain spirals **outward**, but it winds around a fixed
centre: ``S_n - \zeta(s) \approx \frac{n^{1-s}}{1-s}``. The centre is exactly ``\zeta(s)``. Subtracting the
Euler–Maclaurin correction (`partial_sums(s, N; corrected = true)`) turns the outward spiral into one that winds
**inward** onto the centre. At a nontrivial zero that centre is the origin:

```@example ch5
plot_partial_sum_spiral()
```

Zooming in on the origin itself, the curve ``t \mapsto \zeta(\sigma + it)`` passes **through** 0 only for
``\sigma = \tfrac12``, once for every zero. The lines ``\sigma = 0.4`` and ``0.6`` loop around the origin but never
touch it, at least at these heights:

```@example ch5
plot_zeta_near_origin()
```

Putting the two together as an animation, the left panel traces ``\zeta(\tfrac12 + it)`` while the right panel shows
the vector-chain spiral at the current ``t``. Its centre is the orange point, which passes through the origin at each
zero:

```@example ch5
record_zeta_spiral("zeta_spiral.gif"; tmax = 40, frames = 120, framerate = 12)
nothing # hide
```

![animated spiral of partial sums along the critical line](zeta_spiral.gif)

## The Riemann–Siegel formula

Computing ``\zeta(\tfrac12 + it)`` by Euler–Maclaurin costs ``O(t)`` operations. Siegel (1932), working from
Riemann's unpublished notes, found

```math
Z(t) = 2\sum_{n \le N} \frac{\cos(\theta(t) - t\log n)}{\sqrt n}
     + (-1)^{N-1}\left(\frac{t}{2\pi}\right)^{-1/4} C_0(p) + O(t^{-3/4}),
```

with ``N = \lfloor\sqrt{t/2\pi}\rfloor``, ``p = \{\sqrt{t/2\pi}\}`` and
``C_0(p) = \frac{\cos 2\pi(p^2 - p - 1/16)}{\cos 2\pi p}``. Only ``O(\sqrt t)`` terms are needed, and every
large-scale verification of RH is built on this formula.

```@example ch5
riemann_siegel_Z(10_000.0), hardy_Z(10_000.0)
```

## Counting zeros: N(T)

Let ``N(T) = \#\{\rho : 0 < \gamma \le T\}``. The argument principle applied to ``\xi`` around the rectangle
``[-1, 2] \times [0, T]`` gives the **Riemann–von Mangoldt formula**

```math
N(T) = \frac{\theta(T)}{\pi} + 1 + S(T)
     = \frac{T}{2\pi}\log\frac{T}{2\pi e} + \frac78 + S(T) + O(1/T),
```

where ``S(T) = \frac1\pi \arg\zeta(\tfrac12 + iT)`` is defined by continuous variation from ``+\infty + iT``.
``S(T)`` is small (``O(\log T)``) and oscillates wildly. [`argzeta_S`](@ref) computes it by tracking the
argument along the horizontal segment, and [`zero_count`](@ref) rounds the total.

```@example ch5
zero_count(100), zero_count(1000)
```

```@example ch5
plot_zero_counting()
```

## Verifying RH up to height T

If the number of sign changes of ``Z`` on ``(0, T]`` equals ``N(T)``, then **every** zero up to height ``T``
is on the critical line. This is the logic used by Turing, and later by Odlyzko, Gourdon, and Platt–Trudgian
(who reached ``3 \cdot 10^{12}``).

```@example ch5
check_zeros(1000)
```

## Gram points and Gram's law

Gram points ``g_n`` solve ``\theta(g_n) = n\pi``. Since ``Z(t) = \operatorname{Re}\big(e^{i\theta}\zeta\big)`` and
``\zeta(\tfrac12+it) \approx 1 + \dots``, one expects ``(-1)^n Z(g_n) > 0``. This is **Gram's law**. It
holds surprisingly often, but not always. The first failure is at ``n = 126``:

```@example ch5
gram_point(0), gram_point(126)
```

```@example ch5
gram_law_violations(300)
```

## Is the band of fixed width?

There are two different "bands" here, and they behave differently.

**1. The critical strip ``0 < \operatorname{Re} s < 1`` has fixed width 1.** It is not a statement about
zeros. It is simply the region that the two elementary facts leave open: the Euler product rules out
``\operatorname{Re} s > 1``, and the functional equation mirrors that to ``\operatorname{Re} s < 0``. These
boundaries do not depend on the height ``t``.

**2. The band where zeros could still be hiding is not of fixed width.** What has actually been proved depends
on ``t``:

* **Up to ``t = 3\cdot10^{12}`` the band has width 0.** Platt and Trudgian (2021) verified, by the
  sign-change and counting method above ([`check_zeros`](@ref) is a toy version), that every zero there lies
  exactly on the line ([`RH_VERIFIED_HEIGHT`](@ref)).
* **Above that height, only zero-free regions narrow the strip.** For example
  ```math
  \zeta(\sigma+it) \neq 0 \quad\text{for}\quad \sigma \ge 1 - \frac{1}{5.573412\,\log t}\qquad (t \ge 2)
  ```
  (Mossinghoff–Trudgian 2015), and by symmetry also for ``\sigma \le \frac{1}{5.573412\log t}``. The excluded
  slivers have width about ``1/\log t``, which **shrinks to 0** as ``t \to \infty``. So the band where zeros are
  not yet excluded *widens* towards the full width 1 the higher you go. The Vinogradov–Korobov region,
  ``1-\sigma \gg (\log t)^{-2/3}(\log\log t)^{-1/3}``, shrinks more slowly, but it still shrinks.
  [`zero_free_boundary`](@ref) evaluates both.
* **Statistically, the zeros crowd toward the line.** Bohr and Landau (1914) showed that for every
  ``\varepsilon > 0``, all but an infinitesimal proportion of the zeros lie within ``\varepsilon`` of the
  critical line, and zero-density estimates ([Chapter 8](@ref ch-equivalents)) make this quantitative.
  So *almost all* zeros lie in an arbitrarily thin band, but no proof rules out a sparse set of exceptions
  near the edges.
* **RH says the band is a single line, of width 0, at every height.**

```@example ch5
zero_free_boundary(1e12), zero_free_boundary(1e100)
```

```@example ch5
1 - zero_free_boundary(1e100; method = :vinogradov_korobov)
```

```@example ch5
plot_strip_width()
```

## Close pairs (Lehmer's phenomenon)

Occasionally two zeros come extremely close together, and ``Z`` barely crosses the axis between them.
These near-misses are where a counterexample to RH would most plausibly hide. They also drive the
lower bound on the de Bruijn–Newman constant ([Chapter 9](@ref ch-dbn)).

```@example ch5
close_pairs(nontrivial_zeros(3000); k = 3)
```

## Functions in this chapter

[`KNOWN_ZEROS`](@ref), [`riemann_siegel_theta`](@ref), [`hardy_Z`](@ref), [`riemann_siegel_Z`](@ref),
[`gram_point`](@ref), [`gram_points`](@ref), [`gram_law_violations`](@ref), [`argzeta_S`](@ref),
[`riemann_von_mangoldt`](@ref), [`zero_count`](@ref), [`zeros_between`](@ref),
[`nontrivial_zeros`](@ref), [`check_zeros`](@ref), [`close_pairs`](@ref), [`RH_VERIFIED_HEIGHT`](@ref),
[`zero_free_boundary`](@ref), [`partial_sums`](@ref).
