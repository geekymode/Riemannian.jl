# Numerical methods

This page records *how* Riemannian computes things, and how accurately.

## Bernoulli numbers

[`bernoulli`](@ref) uses the recurrence ``\sum_{k=0}^{m}\binom{m+1}{k}B_k = 0`` in exact
`Rational{BigInt}` arithmetic and caches the results. Coefficients such as ``B_{2k}/(2k)!`` are converted
once per floating-point type and precision, and cached as well.

## Complex log-Gamma

[`loggamma_c`](@ref) shifts ``z`` upward with ``\log\Gamma(z) = \log\Gamma(z+m) - \sum_{k<m}\log(z+k)`` until
``|z+m| \ge R``, then applies Stirling's series

```math
\log\Gamma(w) = \left(w - \tfrac12\right)\log w - w + \tfrac12\log 2\pi + \sum_{k=1}^{M} \frac{B_{2k}}{2k(2k-1)w^{2k-1}} .
```

For ``\operatorname{Re} z \ge \tfrac12`` the result is the analytic branch, because each ``\log(z+k)`` has its argument
in ``(-\pi/2, \pi/2)``. This matters for ``\theta(t)``. For ``\operatorname{Re} z < \tfrac12`` the reflection
formula is used, which is correct modulo ``2\pi i``. ``R`` and ``M`` grow with the precision.

## ζ(s)

* **``\operatorname{Re} s \ge 0``:** Euler–Maclaurin ([`zeta_em`](@ref)) with
  ``M = \lfloor p/2\rfloor + 5`` correction terms (``p`` = mantissa bits) and
  ``N = \lceil(|s| + 2M)/\pi\rceil``. Successive correction terms shrink roughly by the factor
  ``|s+2k|^2/(2\pi N)^2 \le 1/4``, so the error is close to the unit roundoff.
* **``\operatorname{Re} s < 0``:** ``\zeta(s) = \chi(s)\zeta(1-s)``, with ``\log\chi`` assembled from
  logarithms. ``\log\sin z`` is evaluated as ``-iz + \log(1-e^{2iz}) + \log(i/2)`` for ``\operatorname{Im} z > 1``
  (and symmetrically for ``\operatorname{Im} z < -1``), so large ``|t|`` does not overflow. Negative integers are
  returned exactly.
* **Cost** is ``O(|s|)``. Against SpecialFunctions.jl the relative error is about ``10^{-16}`` near the real axis and
  about ``10^{-13}`` at ``t = 1000``.

```@example num
using Riemannian
import SpecialFunctions
s = 0.5 + 1000im
abs(zeta(s) - SpecialFunctions.zeta(s)) / abs(zeta(s))
```

## Hardy Z and Riemann–Siegel

[`hardy_Z`](@ref) is computed as ``\operatorname{Re}\left(e^{i\theta(t)}\zeta(\tfrac12+it)\right)``, with ``\zeta`` from
Euler–Maclaurin, so it is accurate but costs ``O(t)``. [`riemann_siegel_Z`](@ref) uses the main sum plus the
``C_0`` correction. It costs ``O(\sqrt t)`` and its error is ``O(t^{-3/4})``.

## Finding and certifying zeros

1. **Scan:** evaluate ``Z`` on a grid with about 4 points per mean gap ``2\pi/\log(t/2\pi)``.
2. **Refine:** for every sign change, run the Illinois (modified regula falsi) method to a relative
   accuracy of ``10^{-13}``.
3. **Certify:** in blocks of about 100 zeros, compare the count with ``N(T_2) - N(T_1)``, where
   ``N(T) = \operatorname{round}\left(\theta(T)/\pi + 1 + S(T)\right)``. ``S(T)`` is obtained by tracking
   ``\arg\zeta(\sigma + iT)`` continuously from ``\sigma = 3`` to ``\tfrac12``, halving the step whenever the
   argument jumps by more than ``\pi/4``.
4. **Repair:** if the counts disagree, for example at a Lehmer pair, rescan that block with 2×, 4×, … finer
   grids.

This is a numerical version of Turing's method. It is *not* interval-certified: no floating-point bounds are
propagated.

## Explicit formulas

``\operatorname{li}(x^\rho)`` is evaluated as ``\operatorname{Ei}(\rho\log x) = -E_1(-\rho\log x) \pm i\pi``, using the
sign of ``\operatorname{Im}``, with SpecialFunctions' complex `expint`. The tail integral in ``J(x)`` is computed
by Simpson's rule after substituting ``u = \log t``.

## de Bruijn–Newman

``\Phi(u)`` is even and analytic in the strip ``|\operatorname{Im} u| < \pi/8``, and it decays doubly
exponentially. For such integrands the trapezoidal rule converges like ``e^{-2\pi d/h}``. With
``h = 1/200`` and 128-bit `BigFloat`s, ``H_t(z)`` is accurate well beyond ``z = 100``.

## Random matrices

GUE matrices are ``H = (A + A^\dagger)/2``, with ``A_{ij}`` standard complex normal. The spectrum fills the
semicircle of radius ``\sqrt{2N}`` and is unfolded with the semicircle's cumulative distribution. Only the central half
of the spectrum is used.
