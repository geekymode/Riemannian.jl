module RiemannianCairoMakieExt

using Riemannian
using CairoMakie
using CairoMakie: Colors
import Riemannian: riemann_theme, plot_prime_counting, plot_prime_gaps, plot_chebyshev,
    plot_mertens, plot_euler_product, plot_zeta_real, plot_analytic_continuation,
    plot_domain_coloring, plot_modulus_surface, plot_critical_line, plot_zeta_spiral,
    plot_zero_counting, plot_explicit_formula, plot_prime_staircase, plot_spacing_distribution,
    plot_pair_correlation, plot_xi, plot_robin, plot_zero_density, plot_heat_flow,
    plot_timeline, gallery

# ── palette & theme ──────────────────────────────────────────────────────────
# Categorical slots in fixed order (validated for colour-vision deficiency on adjacent pairs).
const SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
const INK = "#0b0b0b"         # primary text / reference curves
const INK2 = "#52514e"        # secondary text
const MUTED = "#8a8984"       # axes, guides
const GRID = "#e8e7e3"
const SURFACE = "#fcfcfb"
const SEQ = ["#cde2fb", "#86b6ef", "#3987e5", "#256abf", "#184f95", "#0d366b"]  # sequential blue

function riemann_theme()
    Theme(
        size = (900, 560),
        fontsize = 14,
        backgroundcolor = SURFACE,
        textcolor = INK,
        palette = (color = SERIES,),
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

# Wegert-style phase portrait: hue = arg ζ(s), brightness shows log₂|ζ| contour bands.
function _phase_color(z)
    (isfinite(z) && !iszero(z)) || return Colors.RGB(0.0, 0.0, 0.0)
    h = mod(rad2deg(angle(z)), 360)
    m = log2(abs(z))
    v = 0.72 + 0.28 * (m - floor(m))
    return Colors.RGB(Colors.HSV(h, 0.85, v))
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
                  subtitle = "hue = arg ζ(s); brightness bands = |ζ| doubling", xlabel = "Re s", ylabel = "Im s",
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
        text!(axk, [1.05, -1.05, 0, 0], [0, 0, 1.05, -1.05]; text = ["0", "π", "π/2", "−π/2"],
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
        surface!(ax, xs, ys, Zs; colormap = SEQ, shading = NoShading)
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
            lines!(ax, real.(w), imag.(w); color = ts, colormap = SEQ[2:end], linewidth = 1.5)
            scatter!(ax, [0], [0]; color = SERIES[2], markersize = 10)
            k == 2 && Colorbar(fig[1, 3]; colormap = SEQ[2:end], limits = (0, tmax), label = "t")
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
            scatter!(ax, zs ./ 2, fill(t, length(zs)); color = SERIES[1], markersize = 7)
        end
        γ = filter(<(zmax / 2), nontrivial_zeros(30))
        scatter!(ax, γ, zero(γ); color = SERIES[2], markersize = 11, label = "t = 0: zeros of ζ")
        axislegend(ax; position = :rb)
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
    ]
    paths = String[]
    for (name, f) in plots
        path = joinpath(dir, name * ".png")
        t = @elapsed save(path, f(); px_per_unit)
        @info "saved $path ($(round(t; digits = 1)) s)"
        push!(paths, path)
    end
    return paths
end

end # module
