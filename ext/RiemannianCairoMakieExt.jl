module RiemannianCairoMakieExt

using Riemannian
using CairoMakie
using CairoMakie: Colors
import Riemannian: riemann_theme, plot_prime_counting, plot_prime_gaps, plot_chebyshev,
    plot_mertens, plot_euler_product, plot_zeta_real, plot_analytic_continuation,
    plot_domain_coloring, plot_modulus_surface, plot_critical_line, plot_zeta_spiral,
    plot_zero_counting, plot_explicit_formula, plot_prime_staircase, plot_spacing_distribution,
    plot_pair_correlation, plot_xi, plot_robin, plot_zero_density, plot_heat_flow,
    plot_timeline, plot_critical_strip, plot_strip_schematic, plot_strip_width,
    plot_partial_sum_spiral, plot_zeta_near_origin, record_zeta_spiral, gallery

# ── palette & theme ──────────────────────────────────────────────────────────
# Categorical slots in fixed order (validated for colour-vision deficiency on adjacent pairs).
const SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
const INK = "#0b0b0b"         # primary text / reference curves
const INK2 = "#52514e"        # secondary text
const MUTED = "#8a8984"       # axes, guides
const GRID = "#e8e7e3"
const SURFACE = "#fcfcfb"
const CMAP = :viridis        # every continuous colour scale
const VIRIDIS = Makie.to_colormap(CMAP)

function riemann_theme()
    Theme(
        size = (900, 560),
        fontsize = 14,
        backgroundcolor = SURFACE,
        textcolor = INK,
        palette = (color = SERIES,),
        colormap = CMAP,
        Axis = (
            backgroundcolor = SURFACE,
            xgridcolor = GRID, ygridcolor = GRID,
            topspinevisible = false, rightspinevisible = false,
            leftspinecolor = MUTED, bottomspinecolor = MUTED,
            xtickcolor = MUTED, ytickcolor = MUTED,
            xticklabelcolor = INK2, yticklabelcolor = INK2,
            xlabelcolor = INK2, ylabelcolor = INK2,
            titlealign = :left, titlesize = 17, subtitlecolor = INK2, subtitlesize = 13,
        ),
        Axis3 = (xlabelcolor = INK2, ylabelcolor = INK2, zlabelcolor = INK2, titlesize = 17),
        Lines = (linewidth = 2,),
        Stairs = (linewidth = 2,),
        Scatter = (markersize = 8, strokecolor = SURFACE, strokewidth = 1),
        Legend = (framevisible = false, labelcolor = INK2, patchsize = (22, 12)),
    )
end

_themed(f) = with_theme(f, riemann_theme())

# ── Chapter 1: primes ────────────────────────────────────────────────────────

function plot_prime_counting(; xmax = 100, xbig = 10^6)
    _themed() do
        fig = Figure(size = (1000, 480))
        ax = Axis(fig[1, 1]; title = "Counting primes", subtitle = "π(x) and its classical approximations",
                  xlabel = "x", ylabel = "count")
        tbl = primepi_table(xmax)
        xs = range(2, xmax; length = 600)
        stairs!(ax, 1:xmax, tbl; color = INK, step = :post, label = "π(x)")
        lines!(ax, xs, x -> x / log(x); color = SERIES[1], label = "x / log x")
        lines!(ax, xs, li; color = SERIES[2], label = "li(x)")
        lines!(ax, xs, riemann_R; color = SERIES[3], label = "R(x)")
        axislegend(ax; position = :lt)

        ax2 = Axis(fig[1, 2]; title = "Errors", subtitle = "approximation − π(x), x ≤ $(xbig)",
                   xlabel = "x", ylabel = "error", xscale = log10)
        big = primepi_table(xbig)
        xb = unique(round.(Int, 10 .^ range(1, log10(xbig); length = 400)))
        lines!(ax2, xb, [x / log(x) - big[x] for x in xb]; color = SERIES[1], label = "x / log x")
        lines!(ax2, xb, [li(x) - big[x] for x in xb]; color = SERIES[2], label = "li(x)")
        lines!(ax2, xb, [riemann_R(x) - big[x] for x in xb]; color = SERIES[3], label = "R(x)")
        hlines!(ax2, 0; color = MUTED, linestyle = :dash, linewidth = 1)
        axislegend(ax2; position = :lb)
        fig
    end
end

function plot_prime_gaps(; N = 10^6)
    _themed() do
        ps, gs = prime_gaps(N)
        fig = Figure()
        ax = Axis(fig[1, 1]; title = "Prime gaps", subtitle = "gap to the next prime, p ≤ $N",
                  xlabel = "p", ylabel = "gap", xscale = log10)
        scatter!(ax, ps, gs; color = (SERIES[1], 0.25), markersize = 3, strokewidth = 0,
                 label = "gₙ = pₙ₊₁ − pₙ")
        # record gaps
        rec = Int[]; best = 0
        for (i, g) in enumerate(gs)
            g > best && (best = g; push!(rec, i))
        end
        scatter!(ax, ps[rec], gs[rec]; color = SERIES[2], markersize = 9, label = "record gaps")
        xs = 10 .^ range(0.5, log10(N); length = 200)
        lines!(ax, xs, x -> log(x)^2; color = INK, label = "log² p (Cramér)")
        lines!(ax, xs, log; color = INK2, linestyle = :dash, label = "log p (average gap)")
        axislegend(ax; position = :lt)
        fig
    end
end

function plot_chebyshev(; xmax = 200)
    _themed() do
        fig = Figure()
        ax = Axis(fig[1, 1]; title = "Chebyshev functions",
                  subtitle = "ψ(x) = Σ Λ(n) and ϑ(x) = Σ log p both track x (Prime Number Theorem)",
                  xlabel = "x", ylabel = "")
        xs = 1:xmax
        lines!(ax, [0, xmax], [0, xmax]; color = MUTED, linestyle = :dash, linewidth = 1.5, label = "y = x")
        stairs!(ax, xs, chebyshev_psi.(xs); color = SERIES[1], step = :post, label = "ψ(x)")
        stairs!(ax, xs, chebyshev_theta.(xs); color = SERIES[2], step = :post, label = "ϑ(x)")
        axislegend(ax; position = :lt)
        fig
    end
end

function plot_mertens(; N = 10^5)
    _themed() do
        M = mertens_table(N)
        fig = Figure()
        ax = Axis(fig[1, 1]; title = "Mertens function",
                  subtitle = L"M(x) = \sum_{n\leq x} \mu(n);\quad \mathrm{RH} \Leftrightarrow M(x) = O(x^{1/2+\varepsilon})",
                  xlabel = "x", ylabel = "M(x)")
        xs = 1:N
        band!(ax, xs, -sqrt.(xs), sqrt.(xs); color = (SERIES[1], 0.08))
        lines!(ax, xs, sqrt.(xs); color = SERIES[1], linestyle = :dash, linewidth = 1.5, label = "±√x")
        lines!(ax, xs, -sqrt.(xs); color = SERIES[1], linestyle = :dash, linewidth = 1.5)
        lines!(ax, xs, M; color = INK, linewidth = 1, label = "M(x)")
        axislegend(ax; position = :lb)
        fig
    end
end

# ── Chapter 2: zeta, Euler product, continuation ─────────────────────────────

function plot_euler_product(; s = 2, N = 200)
    _themed() do
        ζ = zeta(s)
        fig = Figure()
        ax = Axis(fig[1, 1]; title = "Sum vs product",
                  subtitle = "error in approximating ζ($s) using integers ≤ N, or primes ≤ N",
                  xlabel = "N", ylabel = "|error|", yscale = log10, xscale = log10)
        Ns = 2:N
        lines!(ax, Ns, [abs(dirichlet_partial(s, n) - ζ) for n in Ns]; color = SERIES[1],
               label = L"\sum_{n\leq N} n^{-s}")
        lines!(ax, Ns, [abs(euler_product_partial(s, n) - ζ) for n in Ns]; color = SERIES[2],
               label = L"\prod_{p\leq N} (1 - p^{-s})^{-1}")
        axislegend(ax; position = :rt)
        fig
    end
end

function plot_zeta_real(; σmin = -12, σmax = 4)
    _themed() do
        fig = Figure(size = (1000, 480))
        ax = Axis(fig[1, 1]; title = "ζ(σ) on the real axis", subtitle = "simple pole at σ = 1",
                  xlabel = "σ", ylabel = "ζ(σ)", limits = (-2, σmax, -4, 6))
        σl = range(-2, 0.97; length = 400); σr = range(1.03, σmax; length = 400)
        lines!(ax, σl, zeta.(σl); color = SERIES[1])
        lines!(ax, σr, zeta.(σr); color = SERIES[1])
        vlines!(ax, 1; color = MUTED, linestyle = :dash, linewidth = 1)
        hlines!(ax, 0; color = MUTED, linewidth = 1)
        pts = [(0.0, -0.5, "ζ(0) = −1/2"), (-1.0, -1 / 12, "ζ(−1) = −1/12"), (2.0, π^2 / 6, "ζ(2) = π²/6")]
        scatter!(ax, first.(pts), getindex.(pts, 2); color = SERIES[2])
        text!(ax, first.(pts), getindex.(pts, 2); text = last.(pts), offset = (6, 8), fontsize = 12,
              color = INK2)

        ax2 = Axis(fig[1, 2]; title = "Trivial zeros", subtitle = "zoom on σ < −1: ζ(−2n) = 0, ζ(−3) = 1/120, ζ(−5) = −1/252, …",
                   xlabel = "σ", ylabel = "ζ(σ)", limits = (σmin, -1, -0.03, 0.03))
        σs = range(σmin, -1; length = 1500)
        hlines!(ax2, 0; color = MUTED, linewidth = 1)
        lines!(ax2, σs, zeta.(σs); color = SERIES[1])
        tz = trivial_zeros(floor(Int, -σmin / 2))
        scatter!(ax2, tz, zero(tz); color = SERIES[2], markersize = 10, label = "−2, −4, −6, …")
        axislegend(ax2; position = :lt)
        fig
    end
end

function plot_analytic_continuation(; σmin = -3, σmax = 3)
    _themed() do
        fig = Figure()
        ax = Axis(fig[1, 1]; title = "Analytic continuation along the real axis",
                  subtitle = "series that converge only on part of the line vs the continued ζ(σ)",
                  xlabel = "σ", ylabel = "value", limits = (σmin, σmax, -4, 6))
        σs = [σ for σ in range(σmin, σmax; length = 700) if abs(σ - 1) > 0.02]
        seg(f) = [abs(σ - 1) < 0.03 ? NaN : f(σ) for σ in σs]
        vlines!(ax, [0, 1]; color = MUTED, linestyle = :dot, linewidth = 1)
        lines!(ax, σs, [real(dirichlet_partial(σ, 10)) for σ in σs]; color = SERIES[1], label = L"\sum_{n\leq 10} n^{-\sigma}")
        lines!(ax, σs, [real(dirichlet_partial(σ, 1000)) for σ in σs]; color = SERIES[3], label = L"\sum_{n\leq 1000} n^{-\sigma}")
        lines!(ax, σs, seg(σ -> real(eta_partial(σ, 1000)) / (1 - 2^(1 - σ))); color = SERIES[2],
               label = L"\eta_{1000}(\sigma) / (1 - 2^{1-\sigma})")
        lines!(ax, σs, seg(zeta); color = INK, linewidth = 2.5, linestyle = :dash, label = "ζ(σ) (continued)")
        text!(ax, [0.02, 1.02], [5.5, 5.5]; text = ["η converges →", "Σ n⁻σ converges →"],
              fontsize = 11, color = INK2)
        axislegend(ax; position = :lb)
        fig
    end
end

# Phase portrait: colour = arg ζ(s) through viridis (jump at arg = ±π, so the branch cuts
# from each zero show as sharp edges), brightness bands mark doublings of |ζ|.
function _phase_color(z)
    (isfinite(z) && !iszero(z)) || return Colors.RGB(0.0, 0.0, 0.0)
    x = (angle(z) + π) / 2π
    c = VIRIDIS[clamp(round(Int, x * (length(VIRIDIS) - 1)) + 1, 1, length(VIRIDIS))]
    m = log2(abs(z))
    v = 0.62 + 0.38 * (m - floor(m))
    return Colors.RGB(v * Colors.red(c), v * Colors.green(c), v * Colors.blue(c))
end

function plot_domain_coloring(; re = (-9, 10), im = (-4, 36), n = 600, nzeros = 20)
    _themed() do
        nx = n
        ny = round(Int, n * (im[2] - im[1]) / (re[2] - re[1]))
        xs = range(re[1], re[2]; length = nx)
        ys = range(im[1], im[2]; length = ny)
        img = Matrix{Colors.RGB{Float64}}(undef, nx, ny)
        Threads.@threads for j in 1:ny
            for i in 1:nx
                s = complex(xs[i], ys[j])
                img[i, j] = _phase_color(s == 1 ? complex(Inf) : zeta(s))
            end
        end
        fig = Figure()
        ax = Axis(fig[1, 1]; title = "Phase portrait of ζ(s)",
                  subtitle = "colour = arg ζ(s) (viridis); brightness bands = |ζ| doubling", xlabel = "Re s", ylabel = "Im s",
                  xticks = -20:4:20, yticks = -40:10:80)
        image!(ax, re[1] .. re[2], im[1] .. im[2], img)
        limits!(ax, re..., im...)
        vlines!(ax, [0, 1]; color = (:white, 0.8), linewidth = 1, linestyle = :dash)
        vlines!(ax, 0.5; color = :white, linewidth = 1.5)
        γ = filter(g -> im[1] ≤ g ≤ im[2], nontrivial_zeros(nzeros))
        γ = vcat(γ, filter(g -> im[1] ≤ g ≤ im[2], -nontrivial_zeros(nzeros)))
        scatter!(ax, fill(0.5, length(γ)), γ; color = :white, strokecolor = INK, strokewidth = 1.5,
                 markersize = 9, label = "nontrivial zeros (Re s = 1/2)")
        tz = filter(x -> re[1] ≤ x, trivial_zeros(20))
        scatter!(ax, tz, zero(tz); color = :white, strokecolor = INK, strokewidth = 1.5, marker = :rect,
                 markersize = 9, label = "trivial zeros")
        scatter!(ax, [1], [0]; color = :black, strokecolor = :white, strokewidth = 1.5, marker = :xcross,
                 markersize = 12, label = "pole s = 1")
        rowsize!(fig.layout, 1, Fixed(820))
        colsize!(fig.layout, 1, Aspect(1, (re[2] - re[1]) / (im[2] - im[1])))
        side = GridLayout(fig[1, 2]; tellheight = false)
        Legend(side[1, 1], ax; framevisible = false)
        axk = Axis(side[2, 1]; width = 140, height = 140, aspect = 1, title = "arg ζ key", titlesize = 12)
        hidedecorations!(axk); hidespines!(axk)
        r = range(-1, 1; length = 160)
        key = [abs(complex(a, b)) ≤ 1 ? _phase_color(complex(a, b) * 1.5) : Colors.RGB(0.988, 0.988, 0.984)
               for a in r, b in r]
        image!(axk, -1 .. 1, -1 .. 1, key)
        text!(axk, [1.05, -1.05, 0, 0], [0, 0, 1.05, -1.05]; text = ["0", "±π", "π/2", "−π/2"],
              align = [(:left, :center), (:right, :center), (:center, :bottom), (:center, :top)],
              fontsize = 11, color = INK2)
        limits!(axk, -1.6, 1.6, -1.6, 1.6)
        resize_to_layout!(fig)
        fig
    end
end

function plot_modulus_surface(; re = (-1, 2), im = (0, 40), n = 220)
    _themed() do
        xs = range(re[1], re[2]; length = n)
        ys = range(im[1], im[2]; length = 2n)
        Zs = [complex(x, y) == 1 ? 3.0 : clamp(log(abs(zeta(complex(x, y)))), -3.0, 3.0) for x in xs, y in ys]
        fig = Figure(size = (1000, 700))
        ax = Axis3(fig[1, 1]; title = "The landscape log|ζ(s)|",
                   xlabel = "Re s", ylabel = "Im s", zlabel = "log|ζ|",
                   azimuth = 1.25π, elevation = 0.18π, aspect = (1, 2.2, 0.7))
        surface!(ax, xs, ys, Zs; colormap = CMAP, shading = NoShading)
        Colorbar(fig[1, 2]; colormap = CMAP, limits = (-3, 3), label = "log|ζ| (clipped)", height = Relative(0.6))
        lines!(ax, fill(0.5, length(ys)), collect(ys), fill(-3.0, length(ys)); color = SERIES[2],
               linewidth = 2)
        text!(ax, Point3f(0.5, im[2], -3.0); text = "critical line", color = SERIES[2], fontsize = 12)
        fig
    end
end

# ── Chapter 3: critical line & zeros ─────────────────────────────────────────

function plot_critical_line(; tmax = 60)
    _themed() do
        ts = range(0, tmax; length = 3000)
        γ = filter(<(tmax), nontrivial_zeros(ceil(Int, tmax / 2) + 5))
        nmax = 0
        while gram_point(nmax + 1) < tmax; nmax += 1; end
        g = [gram_point(n) for n in -1:nmax]
        g = filter(x -> 0 < x < tmax, g)
        fig = Figure(size = (1000, 640))
        ax = Axis(fig[1, 1]; title = "Hardy's Z-function on the critical line",
                  subtitle = "Z(t) = exp(iθ(t)) ζ(1/2+it) is real; each sign change is a zero ρ = 1/2 + iγ",
                  ylabel = "Z(t)")
        hlines!(ax, 0; color = MUTED, linewidth = 1)
        vlines!(ax, g; color = (SERIES[3], 0.5), linewidth = 1, linestyle = :dot)
        lines!(ax, ts, hardy_Z.(ts); color = SERIES[1], label = "Z(t)")
        scatter!(ax, γ, zero(γ); color = SERIES[2], markersize = 10, label = "zeros γₙ")
        scatter!(ax, g, zero(g); color = SERIES[3], marker = :vline, markersize = 14, strokewidth = 0,
                 label = "Gram points")
        axislegend(ax; position = :lt, orientation = :horizontal)
        ax2 = Axis(fig[2, 1]; ylabel = "|ζ(1/2+it)|", xlabel = "t")
        linkxaxes!(ax, ax2)
        lines!(ax2, ts, t -> abs(zeta(complex(0.5, t))); color = INK)
        scatter!(ax2, γ, zero(γ); color = SERIES[2], markersize = 10)
        rowsize!(fig.layout, 2, Relative(0.35))
        fig
    end
end

function plot_zeta_spiral(; tmax = 50)
    _themed() do
        ts = range(0, tmax; length = 4000)
        fig = Figure(size = (1100, 560))
        for (k, σ) in enumerate((0.5, 0.6))
            w = [zeta(complex(σ, t)) for t in ts]
            ax = Axis(fig[1, k]; title = "t ↦ ζ($σ + it), 0 ≤ t ≤ $tmax",
                      subtitle = σ == 0.5 ? "passes through 0 at every zero" : "off the line: misses 0",
                      xlabel = "Re ζ", ylabel = "Im ζ", aspect = DataAspect())
            lines!(ax, real.(w), imag.(w); color = ts, colormap = CMAP, linewidth = 1.5)
            scatter!(ax, [0], [0]; color = SERIES[2], markersize = 10)
            k == 2 && Colorbar(fig[1, 3]; colormap = CMAP, limits = (0, tmax), label = "t")
        end
        fig
    end
end

function plot_zero_counting(; T = 100)
    _themed() do
        γ = filter(<(T), nontrivial_zeros(round(Int, T / 2) + 10))
        ts = range(1, T; length = 2000)
        Nt = [count(<(t), γ) for t in ts]
        smooth = riemann_von_mangoldt.(ts)
        fig = Figure(size = (1000, 640))
        ax = Axis(fig[1, 1]; title = "Counting zeros", subtitle = "N(T) = θ(T)/π + 1 + S(T)",
                  ylabel = "N(t)")
        stairs!(ax, vcat(0, γ, T), vcat(0, 1:length(γ), length(γ)); color = INK, step = :post, label = "N(t)")
        lines!(ax, ts, smooth; color = SERIES[1], label = "θ(t)/π + 1")
        axislegend(ax; position = :lt)
        ax2 = Axis(fig[2, 1]; ylabel = "S(t)", xlabel = "t")
        linkxaxes!(ax, ax2)
        hlines!(ax2, 0; color = MUTED, linewidth = 1)
        lines!(ax2, ts, Nt .- smooth; color = SERIES[2], linewidth = 1.5)
        rowsize!(fig.layout, 2, Relative(0.35))
        fig
    end
end

# ── Chapter 4: explicit formulas ─────────────────────────────────────────────

function plot_explicit_formula(; xmax = 50, ks = (5, 20, 100))
    _themed() do
        γ = nontrivial_zeros(maximum(ks))
        xs = range(1.5, xmax; length = 2500)
        fig = Figure()
        ax = Axis(fig[1, 1]; title = "Primes from zeros: the explicit formula",
                  subtitle = L"\psi(x) = x - \sum_\rho x^\rho/\rho - \log 2\pi - \frac{1}{2}\log(1 - x^{-2}),\ \mathrm{truncated\ to}\ k\ \mathrm{zero\ pairs}",
                  xlabel = "x", ylabel = "ψ(x)")
        for (i, k) in enumerate(ks)
            lines!(ax, xs, x -> psi_explicit(x, γ[1:k]); color = SERIES[i], linewidth = 1.5,
                   label = "k = $k")
        end
        stairs!(ax, 1:xmax, chebyshev_psi.(1:xmax); color = INK, step = :post, label = "ψ(x)")
        axislegend(ax; position = :lt)
        fig
    end
end

function plot_prime_staircase(; xmax = 60, ks = (10, 50, 200), npts = 700)
    _themed() do
        γ = nontrivial_zeros(maximum(ks))
        xs = range(1.6, xmax; length = npts)
        fig = Figure()
        ax = Axis(fig[1, 1]; title = "Riemann's prime staircase",
                  subtitle = L"\pi(x) = \sum_n \frac{\mu(n)}{n} J(x^{1/n}),\quad J(x) = \mathrm{li}(x) - \sum_\rho \mathrm{li}(x^\rho) - \log 2 + \int_x^\infty \frac{dt}{t(t^2-1)\log t}",
                  xlabel = "x", ylabel = "π(x)")
        lines!(ax, xs, riemann_R; color = MUTED, linestyle = :dash, linewidth = 1.5, label = "R(x) (no zeros)")
        for (i, k) in enumerate(ks)
            lines!(ax, xs, x -> primepi_explicit(x, γ[1:k]); color = SERIES[i], linewidth = 1.5,
                   label = "k = $k zeros")
        end
        stairs!(ax, 1:xmax, primepi_table(xmax); color = INK, step = :post, label = "π(x)")
        axislegend(ax; position = :lt)
        fig
    end
end

# ── Chapter 5: statistics ────────────────────────────────────────────────────

function plot_spacing_distribution(; nzeros = 3000, gue = true)
    _themed() do
        δ = normalized_spacings(nontrivial_zeros(nzeros))
        fig = Figure()
        ax = Axis(fig[1, 1]; title = "Spacings between zeros",
                  subtitle = "first $nzeros zeros, unfolded to mean spacing 1",
                  xlabel = "normalised spacing s", ylabel = "density", limits = (0, 3.5, 0, nothing))
        hist!(ax, δ; bins = range(0, 3.5; length = 50), normalization = :pdf,
              color = (SERIES[1], 0.35), strokecolor = SURFACE, strokewidth = 1, label = "ζ zeros")
        if gue
            hist!(ax, gue_spacings(); bins = range(0, 3.5; length = 50), normalization = :pdf,
                  color = (:white, 0), strokecolor = SERIES[7], strokewidth = 1.5,
                  label = "GUE matrices (sampled)")
        end
        ss = range(0, 3.5; length = 300)
        lines!(ax, ss, wigner_gue; color = SERIES[2], label = "GUE surmise")
        lines!(ax, ss, wigner_goe; color = SERIES[3], linestyle = :dash, label = "GOE surmise")
        lines!(ax, ss, poisson_spacing; color = INK2, linestyle = :dot, label = "Poisson (independent)")
        axislegend(ax; position = :rt)
        fig
    end
end

function plot_pair_correlation(; nzeros = 3000, umax = 3)
    _themed() do
        c, d = pair_correlation(nontrivial_zeros(nzeros); umax)
        fig = Figure()
        ax = Axis(fig[1, 1]; title = "Montgomery's pair correlation",
                  subtitle = "density of differences between unfolded zeros",
                  xlabel = "u", ylabel = "density", limits = (0, umax, 0, 1.4))
        barplot!(ax, c, d; color = (SERIES[1], 0.4), gap = 0.1, strokewidth = 0, label = "ζ zeros")
        us = range(0, umax; length = 400)
        lines!(ax, us, montgomery_pair_correlation; color = SERIES[2], label = "1 − (sin πu / πu)²")
        hlines!(ax, 1; color = INK2, linestyle = :dot, linewidth = 1.5, label = "uncorrelated")
        axislegend(ax; position = :rb)
        fig
    end
end

function plot_xi(; tmax = 50)
    _themed() do
        ts = range(0, tmax; length = 2000)
        γ = filter(<(tmax), nontrivial_zeros(ceil(Int, tmax / 2) + 5))
        fig = Figure()
        ax = Axis(fig[1, 1]; title = "Riemann's Ξ function",
                  subtitle = "Ξ(t) = ξ(1/2 + it) is real and even; RH ⇔ all its zeros are real. Scaled by exp(πt/4) to undo its decay.",
                  xlabel = "t", ylabel = L"\Xi(t)\, e^{\pi t/4}")
        hlines!(ax, 0; color = MUTED, linewidth = 1)
        lines!(ax, ts, t -> Xi(t) * exp(π * t / 4); color = SERIES[1])
        scatter!(ax, γ, zero(γ); color = SERIES[2], markersize = 9)
        fig
    end
end

# ── Chapter 6: criteria & bounds ─────────────────────────────────────────────

function plot_robin(; N = 10^5)
    _themed() do
        σ = sigma_table(N)
        ns = 3:N
        r = [σ[n] / (n * log(log(n))) for n in ns]
        eγ = exp(Base.MathConstants.eulergamma)
        bad = robin_violations(N)
        fig = Figure()
        ax = Axis(fig[1, 1]; title = "Robin's criterion",
                  subtitle = "RH ⇔ σ(n) / (n log log n) < e^γ for every n > 5040",
                  xlabel = "n", ylabel = "σ(n) / (n log log n)", xscale = log10, limits = (nothing, (0, 3)))
        scatter!(ax, ns, r; color = (SERIES[1], 0.3), markersize = 2, strokewidth = 0, label = "n")
        hlines!(ax, eγ; color = INK, linestyle = :dash, linewidth = 1.5, label = "e^γ ≈ 1.781")
        vlines!(ax, 5040; color = MUTED, linestyle = :dot, linewidth = 1)
        scatter!(ax, bad, [σ[n] / (n * log(log(n))) for n in bad]; color = SERIES[2], markersize = 8,
                 label = "exceptions (all ≤ 5040)")
        axislegend(ax; position = :rt)
        fig
    end
end

function plot_zero_density()
    _themed() do
        σs = range(0.5, 1; length = 400)
        E = zero_density_exponents.(σs)
        fig = Figure()
        ax = Axis(fig[1, 1]; title = "How many zeros can lie off the line?",
                  subtitle = L"N(\sigma, T) \ll T^{E(\sigma)+\varepsilon}\ \mathrm{counts\ zeros\ with}\ \mathrm{Re}\,\rho \geq \sigma,\ 0<\mathrm{Im}\,\rho \leq T\quad (\mathrm{RH:\ none\ for}\ \sigma > 1/2)",
                  xlabel = "σ", ylabel = "exponent E(σ)", limits = (0.5, 1, 0, 1.1))
        lines!(ax, σs, getfield.(E, :ingham); color = SERIES[1], label = "Ingham (1940)")
        hux = [σ ≥ 0.7 ? e.huxley : NaN for (σ, e) in zip(σs, E)]
        lines!(ax, σs, hux; color = SERIES[2], label = "Huxley (1972)")
        lines!(ax, σs, getfield.(E, :guth_maynard); color = SERIES[3], linewidth = 3, label = "Guth–Maynard (2024)")
        lines!(ax, σs, getfield.(E, :density_hyp); color = INK, linestyle = :dash, label = "Density Hypothesis")
        hlines!(ax, 1; color = MUTED, linestyle = :dot, linewidth = 1, label = "all zeros (~T log T)")
        axislegend(ax; position = :rt)
        fig
    end
end

function plot_heat_flow(; ts = range(-12, 2; length = 36), zmax = 110, step = 0.4)
    _themed() do
        fig = Figure(size = (1000, 600))
        ax = Axis(fig[1, 1]; title = "de Bruijn–Newman heat flow",
                  subtitle = "real zeros of H_t, plotted as z/2 so that t = 0 gives the γₙ. Going down in t, neighbours collide and leave the real axis.\nΛ = last collision time; RH ⇔ Λ ≤ 0; known 0 ≤ Λ ≤ 0.2 (Rodgers–Tao, Platt–Trudgian)",
                  xlabel = "z / 2", ylabel = "t")
        hspan!(ax, 0, 0.2; color = (SERIES[2], 0.15), label = "0 ≤ Λ ≤ 0.2")
        for t in ts
            zs = heat_flow_zeros(t; zmax, step)
            scatter!(ax, zs ./ 2, fill(t, length(zs)); color = fill(t, length(zs)), colormap = CMAP,
                     colorrange = extrema(ts), markersize = 7, strokewidth = 0)
        end
        γ = filter(<(zmax / 2), nontrivial_zeros(30))
        scatter!(ax, γ, zero(γ); color = :white, strokecolor = INK, strokewidth = 1.5, markersize = 11,
                 label = "t = 0: zeros of ζ")
        axislegend(ax; position = :rb)
        Colorbar(fig[1, 2]; colormap = CMAP, limits = extrema(ts), label = "t")
        fig
    end
end

function plot_timeline()
    _themed() do
        bs = breakthroughs()
        cats = [:foundations, :zeros, :computation, :equivalents, :bounds, :analogy]
        markers = [:circle, :diamond, :rect, :utriangle, :dtriangle, :star5]
        n = length(bs)
        fig = Figure(size = (1100, 26n + 140))
        ax = Axis(fig[1, 1]; title = "A timeline of the Riemann Hypothesis",
                  xlabel = "year", limits = (1725, 2035, -n - 1, 0), yticksvisible = false,
                  yticklabelsvisible = false, ygridvisible = false, leftspinevisible = false)
        for (k, c) in enumerate(cats)
            idx = findall(b -> b.category == c, bs)
            scatter!(ax, [bs[i].year for i in idx], -idx; color = SERIES[k], marker = markers[k],
                     markersize = 11, label = string(c))
        end
        for (i, b) in enumerate(bs)
            right = b.year > 1880
            text!(ax, b.year, -i; text = "$(b.year)  $(b.who): $(b.what)", fontsize = 11, color = INK2,
                  align = (right ? :right : :left, :center), offset = (right ? -10 : 10, 0))
        end
        Legend(fig[0, 1], ax; orientation = :horizontal, framevisible = false)
        fig
    end
end

# ── The critical strip: illustrations ────────────────────────────────────────

const _HALO = (glowcolor = (:white, 0.9), glowwidth = 5)

function plot_critical_strip(; tmax = 100, σlims = (-1.5, 2.5), nt = 1400, nσ = 260)
    _themed() do
        ts = range(0, tmax; length = nt)
        σs = range(σlims...; length = nσ)
        L = Matrix{Float64}(undef, nt, nσ)
        Threads.@threads for j in 1:nσ
            for i in 1:nt
                s = complex(σs[j], ts[i])
                L[i, j] = s == 1 ? 3.0 : clamp(log(abs(zeta(s))), -3.0, 3.0)
            end
        end
        γ = filter(<(tmax), nontrivial_zeros(round(Int, tmax / 2) + 10))
        fig = Figure(size = (1400, 560))
        ax = Axis(fig[1, 1]; title = "The critical strip 0 < Re s < 1",
                  subtitle = "colour = log|ζ(σ + it)|. Zeros are the dark wells, every one of them on the line Re s = 1/2",
                  xlabel = "t = Im s", ylabel = "σ = Re s", yticks = [-1, 0, 0.5, 1, 2])
        heatmap!(ax, ts, σs, L; colormap = CMAP, colorrange = (-3, 3))
        band!(ax, [0, tmax], [0, 0], [1, 1]; color = (:white, 0.06))
        hlines!(ax, [0, 1]; color = :white, linewidth = 1.5, linestyle = :dash)
        hlines!(ax, 0.5; color = :white, linewidth = 1)
        scatter!(ax, γ, fill(0.5, length(γ)); color = :white, strokecolor = INK, strokewidth = 1.5,
                 markersize = 10, label = "nontrivial zeros ρ = 1/2 + iγ")
        scatter!(ax, [0], [1]; color = :black, strokecolor = :white, strokewidth = 1.5, marker = :xcross,
                 markersize = 14, label = "pole s = 1")
        x0 = 0.985tmax
        text!(ax, [x0, x0, x0], [1.75, 0.82, -0.9];
              text = ["Re s > 1: Euler product ⇒ no zeros",
                      "Re s = 1: no zeros (⇔ Prime Number Theorem)",
                      "Re s < 0: mirror image (functional equation), only trivial zeros"],
              align = (:right, :center), fontsize = 13, color = INK, _HALO...)
        text!(ax, 0.012tmax, 0.5; text = "critical line\nRe s = 1/2", align = (:left, :center), fontsize = 12,
              color = INK, _HALO...)
        limits!(ax, 0, tmax, σlims...)
        Colorbar(fig[1, 2]; colormap = CMAP, limits = (-3, 3), label = "log|ζ| (clipped)")
        axislegend(ax; position = :lt, backgroundcolor = (:white, 0.85), framevisible = true, framecolor = (:white, 0))
        fig
    end
end

function plot_strip_schematic(; σlims = (-11, 3.5), tmax = 50)
    _themed() do
        γ = filter(<(tmax), nontrivial_zeros(round(Int, tmax / 2) + 10))
        c_left, c_strip, c_right = VIRIDIS[40], VIRIDIS[128], VIRIDIS[215]
        fig = Figure(size = (1000, 820))
        ax = Axis(fig[1, 1]; title = "Where the zeros of ζ live",
                  subtitle = "three regions of the complex plane, shaded with viridis",
                  xlabel = "Re s", ylabel = "Im s", xticks = -10:2:4, yticks = -50:10:50)
        poly!(ax, Rect(σlims[1], -tmax, -σlims[1], 2tmax); color = (c_left, 0.22))
        poly!(ax, Rect(0, -tmax, 1, 2tmax); color = (c_strip, 0.45))
        poly!(ax, Rect(1, -tmax, σlims[2] - 1, 2tmax); color = (c_right, 0.30))
        vlines!(ax, [0, 1]; color = INK2, linewidth = 1, linestyle = :dash)
        vlines!(ax, 0.5; color = INK, linewidth = 1.5)
        hlines!(ax, 0; color = MUTED, linewidth = 1)
        ys = vcat(γ, -γ)
        scatter!(ax, fill(0.5, length(ys)), ys; color = VIRIDIS[1], strokecolor = :white, strokewidth = 1.2,
                 markersize = 9, label = "nontrivial zeros (on Re s = 1/2)")
        tz = filter(>(σlims[1]), trivial_zeros(10))
        scatter!(ax, tz, zero(tz); color = VIRIDIS[200], strokecolor = INK, strokewidth = 1.2, marker = :rect,
                 markersize = 10, label = "trivial zeros −2, −4, …")
        scatter!(ax, [1], [0]; color = :black, strokecolor = :white, strokewidth = 1.5, marker = :xcross,
                 markersize = 13, label = "pole s = 1")
        text!(ax, (σlims[1] - 0) / 2, 0.8tmax; text = "Re s < 0\nζ(s) = χ(s) ζ(1−s)\nonly trivial zeros",
              align = (:center, :center), fontsize = 13, color = INK, _HALO...)
        text!(ax, -0.25, 0.55tmax; text = "critical strip  0 < Re s < 1  →", align = (:right, :center),
              fontsize = 13, color = INK, font = :bold, _HALO...)
        text!(ax, (1 + σlims[2]) / 2, 0.8tmax; text = "Re s > 1\nEuler product\nno zeros",
              align = (:center, :center), fontsize = 13, color = INK, _HALO...)
        text!(ax, 0.5, -0.92tmax; text = "symmetric under s ↦ 1 − s and s ↦ s̄", align = (:center, :center),
              fontsize = 12, color = INK2, _HALO...)
        limits!(ax, σlims..., -tmax, tmax)
        axislegend(ax; position = :lb, backgroundcolor = (:white, 0.9), framevisible = true, framecolor = (:white, 0))
        fig
    end
end

function plot_strip_width(; log10tmax = 40)
    _themed() do
        fig = Figure(size = (1200, 560))
        lt = range(log10(2.0), log10tmax; length = 600)
        Ls = lt .* log(10)
        δ = [Riemannian._zero_free_width(L, :classical) for L in Ls]
        hv = log10(RH_VERIFIED_HEIGHT)
        ax = Axis(fig[1, 1]; title = "Is the band of fixed width?",
                  subtitle = "The strip is always 0 < σ < 1. What we can prove about where zeros may lie depends on t.",
                  xlabel = "log₁₀ t", ylabel = "σ = Re s", limits = (lt[1], log10tmax, -0.02, 1.02))
        band!(ax, [lt[1], log10tmax], [0, 0], [1, 1]; color = (VIRIDIS[128], 0.35),
              label = "not yet excluded (zeros could hide here)")
        band!(ax, lt, 1 .- δ, ones(length(lt)); color = VIRIDIS[220], label = "proven zero-free (and mirror at σ ≈ 0)")
        band!(ax, lt, zeros(length(lt)), δ; color = VIRIDIS[220])
        mask = lt .≤ hv
        band!(ax, lt[mask], zeros(count(mask)), ones(count(mask)); color = (VIRIDIS[30], 0.85),
              label = "t ≤ 3·10¹²: computed, all zeros on the line")
        hlines!(ax, 0.5; color = :white, linewidth = 2)
        text!(ax, hv / 2, 0.5; text = "zeros exactly on σ = 1/2", align = (:center, :bottom), offset = (0, 6),
              color = :white, fontsize = 13)
        text!(ax, (hv + log10tmax) / 2, 0.5; text = "RH: the band is just this line", align = (:center, :bottom),
              offset = (0, 6), fontsize = 13, color = INK, _HALO...)
        δv = Riemannian._zero_free_width(hv * log(10), :classical)
        text!(ax, hv + 0.5, 0.97; text = "↑ proven zero-free sliver: only ≈ $(round(δv; sigdigits = 2)) wide\n   at t = 3·10¹², and shrinking like 1/log t",
              align = (:left, :top), fontsize = 12, color = INK, _HALO...)
        axislegend(ax; position = :rb, backgroundcolor = (:white, 0.9), framevisible = true, framecolor = (:white, 0), labelsize = 12)

        lt2 = range(0.5, 6; length = 600)                     # log10(log10 t) axis: t up to 10^(10^6)
        L2 = (10 .^ lt2) .* log(10)
        ax2 = Axis(fig[1, 2]; title = "Width of the proven zero-free margin",
                   subtitle = "1 − σ₀(t), shrinking to 0: the provable band widens toward width 1",
                   xlabel = "log₁₀ log₁₀ t", ylabel = "1 − σ₀(t)", yscale = log10)
        lines!(ax2, lt2, [Riemannian._zero_free_width(L, :classical) for L in L2]; color = SERIES[1],
               label = "classical: 1/(5.573 log t)  (Mossinghoff–Trudgian)")
        lines!(ax2, lt2, [Riemannian._zero_free_width(L, :vinogradov_korobov) for L in L2]; color = SERIES[2],
               label = "Vinogradov–Korobov (Ford)")
        vlines!(ax2, log10(hv); color = MUTED, linestyle = :dot, linewidth = 1.5)
        text!(ax2, log10(hv), 1e-2; text = " t = 3·10¹²", fontsize = 11, color = INK2)
        axislegend(ax2; position = :lb, labelsize = 12)
        colsize!(fig.layout, 1, Relative(0.55))
        fig
    end
end

# ── Spirals: partial sums and the critical-line curve near 0 ─────────────────

# Path z₀, S₁, S₂, … as points. (Not `vcat(Point2f(…), pts)`: a Point2f is a static vector,
# so vcat would splice in its coordinates instead of prepending a point.)
_path(S, z0 = 0) = pushfirst!(Point2f.(real.(S), imag.(S)), Point2f(real(z0), imag(z0)))

# Colour values for a path of n points. Float, not an integer range: CairoMakie is pathologically
# slow on integer colour vectors.
_ramp(n) = collect(range(0.0, 1.0; length = n))

# A chain of vectors coloured by index (viridis), from z₀ through S₁, S₂, …
function _vector_chain!(ax, S, z0 = 0; linewidth = 1.6)
    pts = _path(S, z0)
    lines!(ax, pts; color = _ramp(length(pts)), colormap = CMAP, linewidth)
end

function plot_partial_sum_spiral(; zero_t = KNOWN_ZEROS[1], other_t = 18.0, N = 400, Nzoom = 4000)
    _themed() do
        fig = Figure(size = (1250, 1050))
        Label(fig[0, 1:2], "Partial sums spiral around ζ(s), and at a zero the spiral closes in on the origin";
              fontsize = 18, font = :bold, halign = :left, tellwidth = false)
        for (row, t) in enumerate((zero_t, other_t))
            s = complex(0.5, t)
            ζs = zeta(s)
            isz = abs(ζs) < 1e-8
            tag = isz ? "s = 1/2 + $(round(t; digits = 4))i (a zero)" : "s = 1/2 + $(t)i (not a zero)"
            mark!(ax, w, lbl) = (scatter!(ax, [real(w)], [imag(w)]; color = isz ? :white : SERIES[2],
                                          strokecolor = INK, strokewidth = 2, markersize = 12);
                                 text!(ax, real(w), imag(w); text = lbl, offset = (10, -14), fontsize = 13,
                                       color = INK, _HALO...))
            origin!(ax) = scatter!(ax, [0], [0]; marker = :cross, color = INK, markersize = 14)

            # 1. Dirichlet series: raw sums spiral outward, corrected sums spiral in
            ax = Axis(fig[row, 1]; title = "Σ n⁻ˢ   $tag", aspect = DataAspect(), xlabel = "Re", ylabel = "Im",
                      subtitle = "raw sums wind outward around ζ(s); corrected sums (orange) wind in")
            _vector_chain!(ax, partial_sums(s, N))
            C = partial_sums(s, N; corrected = true)
            lines!(ax, real.(C), imag.(C); color = SERIES[2], linewidth = 2)
            origin!(ax); mark!(ax, ζs, isz ? "ζ(s) = 0" : "ζ(s)")

            # 2. zoom on the centre: tail n ≥ n₀ ≈ t of the corrected ζ sums. Their error ≈ −(s/12) n^{−s−1}
            #    rotates smoothly (phase −t log n) while shrinking: an inward spiral onto ζ(s).
            Cz = partial_sums(s, Nzoom; corrected = true)
            n0 = max(10, round(Int, t))
            r = 1.15 * maximum(abs.(Cz[n0:end] .- ζs))
            axz = Axis(fig[row, 2]; title = "zoom on the centre (terms $n0…$Nzoom)", aspect = DataAspect(),
                       subtitle = isz ? "corrected sums spiral into the origin" : "corrected sums spiral into ζ(s) ≠ 0",
                       xlabel = "Re", ylabel = "Im",
                       limits = (real(ζs) - r, real(ζs) + r, imag(ζs) - r, imag(ζs) + r))
            tail = Point2f.(real.(Cz[n0:end]), imag.(Cz[n0:end]))
            lines!(axz, tail; color = log.(n0:Nzoom), colormap = CMAP, linewidth = 1.8)   # log n: most terms crowd the centre
            origin!(axz); mark!(axz, ζs, isz ? "0" : "ζ(s)")
        end
        Colorbar(fig[1:2, 3]; colormap = CMAP, limits = (0, 1), label = "progress through the terms (zoom: log scale)",
                 ticks = ([0, 1], ["first", "last"]))
        fig
    end
end

function plot_zeta_near_origin(; tmax = 40, σs = (0.4, 0.5, 0.6), window = 0.5, nt = 8000)
    _themed() do
        ts = range(1, tmax; length = nt)
        γ = filter(<(tmax), nontrivial_zeros(round(Int, tmax / 2) + 5))
        fig = Figure(size = (1200, 600))
        ax = Axis(fig[1, 1]; title = "Zoom on the origin", aspect = DataAspect(),
                  subtitle = "t ↦ ζ(σ + it), 1 ≤ t ≤ $tmax: only σ = 1/2 hits 0, once per zero",
                  xlabel = "Re ζ", ylabel = "Im ζ", limits = (-window, window, -window, window))
        for (k, σ) in enumerate(σs)
            w = [zeta(complex(σ, t)) for t in ts]
            lines!(ax, real.(w), imag.(w); color = σ == 0.5 ? SERIES[1] : (SERIES[k == 1 ? 2 : 3], 0.9),
                   linewidth = σ == 0.5 ? 2 : 1.3, label = "σ = $σ")
        end
        scatter!(ax, [0], [0]; marker = :cross, color = INK, markersize = 16)
        axislegend(ax; position = :lt, backgroundcolor = (:white, 0.9), framevisible = true, framecolor = (:white, 0))

        ax2 = Axis(fig[1, 2]; title = "Distance from 0", subtitle = "|ζ(σ + it)| along each line: only σ = 1/2 reaches 0, at the zeros γₙ (dotted)",
                   xlabel = "t", ylabel = "|ζ(σ + it)|", yscale = log10)
        for (k, σ) in enumerate(σs)
            lines!(ax2, ts, t -> max(abs(zeta(complex(σ, t))), 1e-16);
                   color = σ == 0.5 ? SERIES[1] : SERIES[k == 1 ? 2 : 3], linewidth = σ == 0.5 ? 1.6 : 1.2)
        end
        vlines!(ax2, γ; color = (MUTED, 0.6), linestyle = :dot, linewidth = 1)
        ylims!(ax2, 1e-3, 10)
        fig
    end
end

function record_zeta_spiral(path::AbstractString = "zeta_spiral.gif"; tmax = 50, frames = 150,
                            framerate = 15, N = 300)
    with_theme(riemann_theme()) do
        ts_all = range(0, tmax; length = 4000)
        curve = [zeta(complex(0.5, t)) for t in ts_all]
        tnow = Observable(0.5)
        fig = Figure(size = (1000, 500))
        ax1 = Axis(fig[1, 1]; title = "t ↦ ζ(1/2 + it)", aspect = DataAspect(), xlabel = "Re", ylabel = "Im",
                   limits = (-1.5, 4, -2.5, 2.5))
        upto = @lift findlast(<=($tnow), ts_all)
        trace = @lift Point2f.(real.(curve[1:$upto]), imag.(curve[1:$upto]))
        lines!(ax1, trace; color = SERIES[1], linewidth = 1.5)
        scatter!(ax1, [0], [0]; marker = :cross, color = INK, markersize = 14)
        head = @lift Point2f(reim(zeta(complex(0.5, $tnow)))...)
        scatter!(ax1, head; color = SERIES[2], markersize = 12)
        label = @lift "t = $(round($tnow; digits = 2))"
        text!(ax1, -1.4, 2.3; text = label, fontsize = 14, color = INK)

        ax2 = Axis(fig[1, 2]; title = "partial sums of Σ n⁻ˢ (n ≤ $N)", aspect = DataAspect(),
                   subtitle = "the vector chain spirals around ζ(1/2 + it): through 0 exactly at a zero",
                   xlabel = "Re", ylabel = "Im", limits = (-1.5, 4, -2.5, 2.5))
        chain = @lift _path(partial_sums(complex(0.5, $tnow), N))
        lines!(ax2, chain; color = _ramp(N + 1), colormap = CMAP, linewidth = 1.4)
        scatter!(ax2, [0], [0]; marker = :cross, color = INK, markersize = 14)
        scatter!(ax2, head; color = SERIES[2], strokecolor = :white, strokewidth = 1.5, markersize = 12)
        record(fig, path, range(0.5, tmax; length = frames); framerate) do t
            tnow[] = t
        end
        path
    end
end

# ── gallery ──────────────────────────────────────────────────────────────────

function gallery(dir::AbstractString = "figures"; px_per_unit = 2, quick = false)
    mkpath(dir)
    plots = [
        "01_prime_counting" => () -> plot_prime_counting(),
        "02_prime_gaps" => () -> plot_prime_gaps(),
        "03_chebyshev" => () -> plot_chebyshev(),
        "04_mertens" => () -> plot_mertens(),
        "05_euler_product" => () -> plot_euler_product(),
        "06_zeta_real" => () -> plot_zeta_real(),
        "07_analytic_continuation" => () -> plot_analytic_continuation(),
        "08_domain_coloring" => () -> plot_domain_coloring(; n = quick ? 250 : 600),
        "09_modulus_surface" => () -> plot_modulus_surface(; n = quick ? 80 : 220),
        "10_critical_line" => () -> plot_critical_line(),
        "11_zeta_spiral" => () -> plot_zeta_spiral(),
        "12_zero_counting" => () -> plot_zero_counting(),
        "13_xi" => () -> plot_xi(),
        "14_explicit_formula" => () -> plot_explicit_formula(),
        "15_prime_staircase" => () -> plot_prime_staircase(; npts = quick ? 200 : 700),
        "16_spacing_distribution" => () -> plot_spacing_distribution(; nzeros = quick ? 1000 : 3000),
        "17_pair_correlation" => () -> plot_pair_correlation(; nzeros = quick ? 1000 : 3000),
        "18_robin" => () -> plot_robin(),
        "19_zero_density" => () -> plot_zero_density(),
        "20_heat_flow" => () -> plot_heat_flow(; ts = quick ? range(-12, 2; length = 8) : range(-12, 2; length = 36)),
        "21_timeline" => () -> plot_timeline(),
        "22_critical_strip" => () -> plot_critical_strip(; nt = quick ? 500 : 1400, nσ = quick ? 100 : 260),
        "23_strip_schematic" => () -> plot_strip_schematic(),
        "24_strip_width" => () -> plot_strip_width(),
        "25_partial_sum_spiral" => () -> plot_partial_sum_spiral(),
        "26_zeta_near_origin" => () -> plot_zeta_near_origin(; nt = quick ? 3000 : 8000),
    ]
    paths = String[]
    gif = joinpath(dir, "zeta_spiral.gif")
    t = @elapsed record_zeta_spiral(gif; tmax = 40, frames = quick ? 40 : 120, framerate = 12)
    @info "saved $gif ($(round(t; digits = 1)) s)"
    push!(paths, gif)
    for (name, f) in plots
        path = joinpath(dir, name * ".png")
        t = @elapsed save(path, f(); px_per_unit)
        @info "saved $path ($(round(t; digits = 1)) s)"
        push!(paths, path)
    end
    return paths
end

end # module
