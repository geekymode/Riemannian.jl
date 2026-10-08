# ─────────────────────────────────────────────────────────────────────────────
# Special-function helpers: Bernoulli numbers, complex log-gamma, stable log-sin.
# Everything is generic in the float type so that BigFloat works too.
# ─────────────────────────────────────────────────────────────────────────────

const _BERNOULLI = Rational{BigInt}[1 // 1]          # B₀
const _CACHE_LOCK = ReentrantLock()     # caches are filled lazily, possibly from several threads

"""
    bernoulli(n) -> Rational{BigInt}

Bernoulli number Bₙ (convention B₁ = −1/2). Exact, cached.
B₂ = 1/6, B₄ = −1/30, … — they give ζ(2k) and ζ(1−2k) exactly.
"""
function bernoulli(n::Integer)
    n ≥ 0 || throw(DomainError(n))
    return lock(_CACHE_LOCK) do
        _extend_bernoulli!(n)
        _BERNOULLI[n+1]
    end
end

function _extend_bernoulli!(n)
    while length(_BERNOULLI) ≤ n
        m = length(_BERNOULLI)                       # computing B_m
        if m > 1 && isodd(m)
            push!(_BERNOULLI, 0 // 1)
            continue
        end
        s = zero(Rational{BigInt})
        for k in 0:m-1
            s += binomial(big(m + 1), k) * _BERNOULLI[k+1]
        end
        push!(_BERNOULLI, -s / (m + 1))
    end
end

const _COEFF_CACHE = Dict{Any,Any}()

# B_{2k}/(2k)!  for k = 1..M, converted once per float type/precision.
function _em_coeffs(::Type{T}, M::Int) where {T}
    key = (:em, T, precision(T), M)
    lock(_CACHE_LOCK) do
        get!(_COEFF_CACHE, key) do
            T[T(bernoulli(2k) / factorial(big(2k))) for k in 1:M]
        end
    end::Vector{T}
end

# B_{2k}/(2k(2k−1)) for Stirling's series.
function _stirling_coeffs(::Type{T}, M::Int) where {T}
    key = (:stirling, T, precision(T), M)
    lock(_CACHE_LOCK) do
        get!(_COEFF_CACHE, key) do
            T[T(bernoulli(2k) / (2k * (2k - 1))) for k in 1:M]
        end
    end::Vector{T}
end

"""
    _logsin(z)

log(sin z) modulo 2πi, stable for large |Im z| where sin overflows.
"""
function _logsin(z::Complex{T}) where {T<:AbstractFloat}
    y = imag(z)
    if y > 1
        return -im * z + log1p(-exp(2im * z)) + log(Complex{T}(0, 1) / 2)
    elseif y < -1
        return im * z + log1p(-exp(-2im * z)) + log(Complex{T}(0, -1) / 2)
    else
        return log(sin(z))
    end
end

"""
    loggamma_c(z)

Complex log-Gamma. For Re z ≥ 1/2 this is the analytic branch (continuous in z),
which matters for the Riemann–Siegel theta function; for Re z < 1/2 the reflection
formula is used and the result is only correct modulo 2πi (fine for exp).
"""
function loggamma_c(z::Complex{T}) where {T<:AbstractFloat}
    if real(z) < 0.5
        return log(T(π)) - _logsin(T(π) * z) - loggamma_c(1 - z)
    end
    bits = precision(T)
    R = T(max(17, bits ÷ 3))
    M = max(10, bits ÷ 6)
    acc = zero(Complex{T})
    w = z
    while abs(w) < R
        acc += log(w)
        w += 1
    end
    c = _stirling_coeffs(T, M)
    w⁻¹ = inv(w)
    w⁻² = w⁻¹ * w⁻¹
    s = zero(Complex{T})
    p = w⁻¹
    for k in 1:M
        s += c[k] * p
        p *= w⁻²
    end
    return (w - T(0.5)) * log(w) - w + log(2 * T(π)) / 2 + s - acc
end
loggamma_c(z::Number) = loggamma_c(complex(float(z)))

# Principal branch of Lambert W for x ≥ 0 (used for Gram-point initial guesses).
function _lambertw(x::Real)
    w = x < 3 ? log1p(x) : log(x) - log(log(x))
    for _ in 1:50
        ew = exp(w)
        f = w * ew - x
        wn = w - f / (ew * (w + 1) - (w + 2) * f / (2w + 2))   # Halley
        abs(wn - w) ≤ 1e-15 * (1 + abs(wn)) && return wn
        w = wn
    end
    return w
end
