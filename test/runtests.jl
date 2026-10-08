using Riemannian
using Test

@testset "Riemannian" begin
    @testset "primes & arithmetic" begin
        @test sieve(30) == [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]
        @test primepi(10^6) == 78498
        @test nthprime(10_000) == 104729
        @test isprime(2^61 - 1) && !isprime(2^61 + 1)
        @test all(isprime(p) == (p in sieve(2000)) for p in 1:2000)
        @test factorize(360) == [2 => 3, 3 => 2, 5 => 1]
        @test mobius_table(10) == [mobius(n) for n in 1:10] == [1, -1, -1, 0, -1, 1, -1, 0, 0, 1]
        @test mertens(1000) == 2
        @test totient(36) == 12
        @test divisor_sigma(12) == 28
        @test sigma_table(100) == [divisor_sigma(n) for n in 1:100]
        @test chebyshev_psi(10) ≈ log(2520)                 # lcm(1..10)
        @test riemann_R(1e6) ≈ 78527.3994 atol = 1e-3
        @test li(1e6) ≈ 78627.5491 atol = 1e-3
    end

    @testset "exact values" begin
        @test bernoulli(1) == -1 // 2
        @test bernoulli(12) == -691 // 2730
        @test zeta_even_exact(1) == 1 // 6
        @test zeta_even_exact(2) == 1 // 90
        @test zeta_negint_exact(0) == -1 // 2
        @test zeta_negint_exact(1) == -1 // 12
        @test zeta_negint_exact(4) == 0
    end

    @testset "zeta" begin
        @test zeta(2) ≈ π^2 / 6 rtol = 1e-14
        @test zeta(-1) ≈ -1 / 12
        @test zeta(-4.0) == 0
        @test zeta(0.5) ≈ -1.4603545088095868 rtol = 1e-13
        @test abs(zeta(complex(0.5, KNOWN_ZEROS[1]))) < 1e-12
        # functional equation
        for s in (0.3 + 7im, -2.5 + 20im, 0.8 - 3im)
            @test zeta(s) ≈ chi_factor(s) * zeta(1 - s) rtol = 1e-11
            @test completed_zeta(s) ≈ completed_zeta(1 - s) rtol = 1e-11
        end
        # three independent continuations agree in the strip
        s = 0.3 + 5im
        @test zeta_borwein(s) ≈ zeta(s) rtol = 1e-12
        @test zeta_em(s; N = 50, M = 20) ≈ zeta(s) rtol = 1e-12
        @test dirichlet_eta(1) ≈ log(2)
        @test abs(euler_product_partial(3, 10^4) - zeta(3)) < 1e-8
        # BigFloat
        setprecision(BigFloat, 256) do
            @test abs(zeta(big(2)) - big(π)^2 / 6) < big(10.0)^-70
        end
        @test abs(imag(xi(0.5 + 10im))) < 1e-14
    end

    @testset "zeros" begin
        z = nontrivial_zeros(10)
        @test maximum(abs.(z .- collect(KNOWN_ZEROS))) < 1e-10
        @test trivial_zeros(3) == [-2, -4, -6]
        @test hardy_Z(KNOWN_ZEROS[1]) ≈ 0 atol = 1e-12
        @test riemann_siegel_Z(1000.0) ≈ hardy_Z(1000.0) atol = 1e-3
        @test gram_point(0) ≈ 17.8455995404 atol = 1e-8
        @test first(gram_law_violations(130)) == 126
        @test zero_count(100) == 29
        @test zero_count(1000) == 649
        @test check_zeros(300).all_on_line
        @test length(nontrivial_zeros(1000)) == 1000
        @test issorted(nontrivial_zeros(1000))
        @test zero_free_boundary(1e12) ≈ 1 - 1 / (5.573412 * log(1e12))
        @test zero_free_boundary(10.0^100) > zero_free_boundary(1e12)          # margin shrinks
        @test all(t -> abs(zeta(complex(zero_free_boundary(t), t))) > 0.01, 10.0:10.0:1000.0)
        @test RH_VERIFIED_HEIGHT == 3e12
    end

    @testset "explicit formulas" begin
        γ = nontrivial_zeros(300)
        @test abs(psi_explicit(100.5, γ) - chebyshev_psi(100.5)) < 1
        @test abs(primepi_explicit(100.5, γ) - primepi(100)) < 0.5
        @test riemann_J(100) ≈ 25 + 4 / 2 + 2 / 3 + 2 / 4 + 1 / 5 + 1 / 6
    end

    @testset "statistics" begin
        δ = normalized_spacings(nontrivial_zeros(1000))
        @test 0.95 < sum(δ) / length(δ) < 1.05
        @test montgomery_pair_correlation(1.0) ≈ 1
        c, d = pair_correlation(nontrivial_zeros(2000))
        @test d[1] < 0.2 && 0.8 < sum(d[end-9:end]) / 10 < 1.2   # repulsion at 0, → 1 far out
        g = gue_spacings(200; samples = 5)
        @test 0.9 < sum(g) / length(g) < 1.1
    end

    @testset "criteria" begin
        @test maximum(robin_violations(20_000)) == 5040
        @test isempty(lagarias_violations(5000))
        @test 0 < li_coefficient(1, nontrivial_zeros(1000)) < LI_LAMBDA1
        @test mertens_ratio_max(10^4).ratio < 1
        @test zero_density_exponents(0.75).guth_maynard < zero_density_exponents(0.75).ingham
    end

    @testset "de Bruijn–Newman" begin
        @test debruijn_H(0, 10.0) ≈ Xi(5.0) / 8 rtol = 1e-10
        z = heat_flow_zeros(0.0; zmax = 45)
        @test z ./ 2 ≈ collect(KNOWN_ZEROS[1:2]) atol = 1e-8
    end

    @testset "narrative" begin
        @test length(breakthroughs()) > 20
        @test occursin("−2n", sprint(io -> explain(:trivial_zeros; io)))
    end
end
