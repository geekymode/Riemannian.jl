# ─────────────────────────────────────────────────────────────────────────────
# Plotting API. Implemented in ext/RiemannianCairoMakieExt.jl — run `using CairoMakie`.
# Every function returns a `Figure`; save with `save("file.png", fig)`.
# ─────────────────────────────────────────────────────────────────────────────

const _PLOT_FUNCTIONS = (
    riemann_theme = "Makie theme used by all plots (apply with `set_theme!(riemann_theme())`).",
    plot_prime_counting = "π(x) staircase vs x/log x, li(x), R(x); and their errors at larger x. `(; xmax = 100, xbig = 10^6)`",
    plot_prime_gaps = "Prime gaps vs p with the log²p (Cramér) envelope. `(; N = 10^6)`",
    plot_chebyshev = "ψ(x) and ϑ(x) against x. `(; xmax = 200)`",
    plot_mertens = "Mertens function M(x) with ±√x. `(; N = 10^5)`",
    plot_euler_product = "Convergence of Dirichlet partial sums and Euler partial products to ζ(s). `(; s = 2, N = 200)`",
    plot_zeta_real = "ζ(σ) on the real axis: pole at 1, ζ(0) = −1/2, trivial zeros. `(; σmin = -14, σmax = 4)`",
    plot_analytic_continuation = "Dirichlet series vs η-series vs continued ζ on the real axis. `(; σmin = -3, σmax = 3)`",
    plot_domain_coloring = "Phase portrait of ζ on a rectangle, marking critical strip, zeros and pole. `(; re = (-12, 12), im = (-5, 45), n = 600)`",
    plot_modulus_surface = "3-D landscape of log|ζ(s)| over part of the critical strip. `(; re = (-1, 2), im = (0, 40))`",
    plot_critical_line = "Hardy Z(t) and |ζ(1/2+it)| with zeros and Gram points. `(; tmax = 60)`",
    plot_zeta_spiral = "Image of the critical line t ↦ ζ(1/2+it) (passes through 0 at zeros) vs σ = 0.6. `(; tmax = 50)`",
    plot_zero_counting = "N(T) staircase vs Riemann–von Mangoldt formula, plus S(T). `(; T = 100)`",
    plot_explicit_formula = "ψ(x) and explicit-formula approximations with k zeros. `(; xmax = 50, ks = (5, 20, 100))`",
    plot_prime_staircase = "π(x) reconstructed from zeros via Riemann's formula. `(; xmax = 60, ks = (10, 50, 200))`",
    plot_spacing_distribution = "Normalised zero spacings vs GUE, GOE, Poisson. `(; nzeros = 3000)`",
    plot_pair_correlation = "Pair correlation of zeros vs Montgomery's 1 − (sin πu/πu)². `(; nzeros = 3000)`",
    plot_xi = "Riemann Ξ(t) = ξ(1/2+it), real with real zeros. `(; tmax = 50)`",
    plot_robin = "Robin's ratio σ(n)/(n log log n) against e^γ. `(; N = 10^5)`",
    plot_zero_density = "Zero-density exponents: Ingham, Huxley, Guth–Maynard, density hypothesis.",
    plot_heat_flow = "de Bruijn–Newman: real zeros of H_t moving with t. `(; ts = range(-12, 2; length = 36), zmax = 110)`",
    plot_timeline = "Timeline of milestones from `breakthroughs()`.",
    gallery = "Render every plot into a directory: `gallery(\"figures\")`.",
)

for (name, doc) in pairs(_PLOT_FUNCTIONS)
    @eval begin
        @doc $("    $(name)(...)\n\n$(doc)\n\nRequires `using CairoMakie`.") function $name end
    end
end
