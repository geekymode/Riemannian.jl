using Documenter
using Riemannian
using CairoMakie

CairoMakie.activate!(type = "png", px_per_unit = 1.5)

DocMeta.setdocmeta!(Riemannian, :DocTestSetup, :(using Riemannian); recursive = true)

const EXT = Base.get_extension(Riemannian, :RiemannianCairoMakieExt)

# Source links need a git commit; drop them for local builds of an uncommitted tree.
const HAS_COMMIT = success(pipeline(`git -C $(dirname(@__DIR__)) rev-parse HEAD`; stdout = devnull, stderr = devnull))
const REPO_KW = HAS_COMMIT ? (; repo = Remotes.GitHub("geekymode", "Riemannian.jl")) : (; remotes = nothing)

makedocs(;
    REPO_KW...,
    sitename = "Riemannian.jl",
    authors = "rethna",
    modules = [Riemannian, EXT],
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://geekymode.github.io/Riemannian.jl",
        edit_link = "main",
        assets = ["assets/custom.css"],
        size_threshold = 2^21,
        size_threshold_warn = 2^20,
        example_size_threshold = 0,          # write every figure to its own file
        mathengine = Documenter.KaTeX(),
    ),
    pages = [
        "Home" => "index.md",
        "Getting started" => "getting_started.md",
        "The story" => [
            "1 · Primes" => "chapters/01_primes.md",
            "2 · Euler product & ζ at integers" => "chapters/02_euler.md",
            "3 · Analytic continuation" => "chapters/03_continuation.md",
            "4 · Trivial zeros" => "chapters/04_trivial_zeros.md",
            "5 · Nontrivial zeros" => "chapters/05_nontrivial_zeros.md",
            "6 · Explicit formulas" => "chapters/06_explicit_formulas.md",
            "7 · Random matrices" => "chapters/07_random_matrices.md",
            "8 · Equivalent criteria & bounds" => "chapters/08_equivalents.md",
            "9 · de Bruijn–Newman" => "chapters/09_debruijn_newman.md",
            "10 · History & breakthroughs" => "chapters/10_history.md",
        ],
        "Gallery" => "gallery.md",
        "Numerical methods" => "numerics.md",
        "API reference" => "api.md",
        "References" => "references.md",
    ],
    checkdocs = :exports,
    warnonly = [:missing_docs],
    doctest = true,
)

deploydocs(;
    repo = "github.com/geekymode/Riemannian.jl",
    devbranch = "main",
    push_preview = true,
)
