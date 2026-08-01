-- [Special Functions]
-- Gamma(x) and Log_gamma(x): the gamma function (generalizes factorial
-- to real arguments -- Gamma(n+1) = n! for every nonnegative integer n)
-- and its natural logarithm -- the two foundational special functions
-- almost everything else in this project builds on. Built on top of
-- them, further down this file: Beta(a,b) and Log_beta(a,b)
-- (Gamma(a)*Gamma(b)/Gamma(a+b), for a,b > 0); Gamma_p(a,x) and
-- Gamma_q(a,x), the regularized lower/upper incomplete gamma functions;
-- and Erf(x)/Erfc(x), the (complementary) error function, which turns
-- out to just BE a regularized incomplete gamma function with a=1/2.
-- Bessel_j(n,x) (the Bessel function of the first kind) is independent
-- of all of the above -- its own algorithm, its own header, further
-- down.
--
-- Claude: computed via the Lanczos approximation (g=7, n=9), the same
-- well-known, validated coefficient set used by Numerical Recipes and
-- Boost.Math -- copied here, not derived, the same discipline as the
-- QK21 quadrature tables and Brent-Dekker's algorithm. Verified directly
-- against Python's own math.gamma/math.lgamma (an independent, trusted
-- reference) across a wide spread of positive, negative-non-integer, and
-- large values -- max relative error ~1.6e-13, at the limit of double
-- precision itself.
--
-- Log_gamma is computed FIRST, entirely in log-space, and Gamma is built
-- on top of it (Gamma(x) = sign(x) * exp(Log_gamma(x))) -- deliberately,
-- not the other way around. Computing Gamma directly via the textbook
-- Lanczos formula (t^(x-0.5) * exp(-t), multiplied as separate factors)
-- overflows for x as small as ~150, EVEN THOUGH the true Gamma(150) is
-- still comfortably representable as a double -- the huge and tiny
-- intermediate factors partially cancel in the true product, but
-- computing them as separate numbers before multiplying loses that
-- cancellation and overflows early, the same class of hazard as
-- stable_tanh's own fix elsewhere in this codebase. Caught directly by
-- testing, not assumed.
local Core = require("GABAS-SPECIAL-FUNCTIONS.modules.core")

local LANCZOS_G = 7
local LANCZOS_C = {
    0.99999999999980993,
    676.5203681218851,
    -1259.1392167224028,
    771.32342877765313,
    -176.61502916214059,
    12.507343278686905,
    -0.13857109526572012,
    9.9843695780195716e-6,
    1.5056327351493116e-7,
}
local HALF_LOG_TWO_PI = 0.5 * math.log(2 * math.pi)

local function is_nonpositive_integer(x)
    return x <= 0 and x == math.floor(x)
end

-- Claude: valid only for x >= 0.5 -- Gamma is strictly positive there,
-- so "log of Gamma" needs no absolute-value handling in this branch.
-- Log_gamma below routes anything smaller through the reflection formula
-- instead of calling this directly.
local function log_gamma_positive_branch(x)
    local z = x - 1
    local a = LANCZOS_C[1]
    for i = 1, LANCZOS_G + 1 do
        a = a + LANCZOS_C[i + 1] / (z + i)
    end
    local t = z + LANCZOS_G + 0.5
    return HALF_LOG_TWO_PI + (z + 0.5) * math.log(t) - t + math.log(a)
end

-- Log_gamma(x) -> log(|Gamma(x)|)
--
-- The natural logarithm of the ABSOLUTE VALUE of Gamma(x) -- the
-- standard convention (matching, e.g., C's lgamma, Python's
-- math.lgamma): Gamma(x) is negative on some intervals of the negative
-- reals (e.g. (-2,-1)), where log(Gamma(x)) itself would not be a real
-- number, so every reference implementation returns log|Gamma(x)|
-- instead. Use Gamma(x) directly when you need the sign too and
-- overflow isn't a concern; use Log_gamma when x is large enough that
-- Gamma(x) itself would overflow a double (Gamma(172) already does) but
-- only a log-scale result is actually needed (e.g. ratios of gamma
-- functions, as Beta and the incomplete gamma/beta functions will need).
--
-- Rejects the nonpositive integers (0, -1, -2, ...) -- genuine poles of
-- Gamma, not merely large values -- with a clear error rather than a
-- silent Inf or NaN.
local function Log_gamma(x)
    Core.assert_finite_number(x, "Log_gamma: x")
    assert(not is_nonpositive_integer(x),
        "Log_gamma: x must not be a nonpositive integer (Gamma has a pole there).")
    if x >= 0.5 then
        return log_gamma_positive_branch(x)
    end
    -- Reflection formula: |Gamma(x)| = pi / (|sin(pi*x)| * Gamma(1-x)),
    -- and 1-x > 0.5 here, so the positive branch above applies to it.
    return math.log(math.pi) - math.log(math.abs(math.sin(math.pi * x))) - log_gamma_positive_branch(1 - x)
end

-- Gamma(x) -> the gamma function at x, generalizing factorial to real
-- arguments: Gamma(n+1) = n! for every nonnegative integer n.
--
-- Claude: the sign for x < 0 comes from the same reflection formula
-- Log_gamma uses -- Gamma(1-x) is always positive there (since
-- 1-x > 1), and pi is positive, so the sign of Gamma(x) is exactly the
-- sign of 1/sin(pi*x), i.e. the sign of sin(pi*x) itself.
--
-- For a large enough x (Gamma(172) and beyond), the true value exceeds
-- what a double can represent -- this returns +-math.huge in that case,
-- the same honest IEEE-754 overflow behavior math.exp already has
-- elsewhere in Lua, not a silently wrong finite number. Use Log_gamma
-- instead when x is that large and only a log-scale result is needed.
local function Gamma(x)
    Core.assert_finite_number(x, "Gamma: x")
    assert(not is_nonpositive_integer(x),
        "Gamma: x must not be a nonpositive integer (Gamma has a pole there).")
    local sign = 1
    if x < 0 then
        sign = math.sin(math.pi * x) > 0 and 1 or -1
    end
    return sign * math.exp(Log_gamma(x))
end

-- Log_beta(a, b) -> log(Beta(a,b))
--
-- Beta(a,b) = Gamma(a)*Gamma(b)/Gamma(a+b), for a > 0 and b > 0 -- on
-- this domain Beta is always strictly positive, so (unlike Log_gamma)
-- no absolute-value/sign subtlety is needed here.
--
-- Claude: computed as Log_gamma(a) + Log_gamma(b) - Log_gamma(a+b) --
-- entirely in log-space, never forming Gamma(a), Gamma(b), or
-- Gamma(a+b) as actual numbers. This is not just tidiness: Beta(100,100)
-- is an ordinary, well-behaved positive number (~2.2e-61), but
-- Gamma(100+100) = Gamma(200) already overflows a double on its own
-- (Gamma(172) is the overflow threshold) -- computing the naive
-- Gamma(a)*Gamma(b)/Gamma(a+b) formula directly produces
-- finite*finite/inf = 0, silently wrong by 61 orders of magnitude, not
-- an honest failure. Verified directly: the naive formula really does
-- overflow at a=b=100 (confirmed against Python before writing this),
-- while this log-space formula stays accurate.
local function Log_beta(a, b)
    Core.assert_finite_number(a, "Log_beta: a")
    Core.assert_finite_number(b, "Log_beta: b")
    assert(a > 0, "Log_beta: a must be positive.")
    assert(b > 0, "Log_beta: b must be positive.")
    return Log_gamma(a) + Log_gamma(b) - Log_gamma(a + b)
end

-- Beta(a, b) -> Gamma(a)*Gamma(b)/Gamma(a+b), for a > 0 and b > 0
--
-- The classical Euler integral of the first kind,
-- Beta(a,b) = integral from 0 to 1 of t^(a-1)*(1-t)^(b-1) dt -- verified
-- directly against that defining integral (independent numerical
-- quadrature), not just against its own Gamma-ratio formula.
--
-- Built on Log_beta (exp(Log_beta(a,b))) rather than computed via
-- Gamma(a)*Gamma(b)/Gamma(a+b) directly -- see Log_beta's own header for
-- exactly why that distinction matters here, not just for large a,b.
local function Beta(a, b)
    Core.assert_finite_number(a, "Beta: a")
    Core.assert_finite_number(b, "Beta: b")
    assert(a > 0, "Beta: a must be positive.")
    assert(b > 0, "Beta: b must be positive.")
    return math.exp(Log_beta(a, b))
end

-- Claude: the regularized incomplete gamma functions -- P(a,x) below and
-- Q(a,x) = 1-P(a,x) -- have no single formula that stays numerically
-- accurate everywhere, so (following the classical Numerical Recipes
-- algorithm, "gser"/"gcf"/"gammp"/"gammq") this uses TWO different
-- representations and picks whichever one is actually converging fast
-- at the given (a,x):
--
--   a power series for P(a,x), accurate when x < a+1 (the series' terms
--   shrink quickly there);
--
--   a continued fraction (evaluated via the modified Lentz method, which
--   avoids ever dividing by an intermediate zero) for Q(a,x), accurate
--   when x >= a+1 (where the series above would need very many terms
--   and accumulate rounding error before converging).
--
-- Both are scaled by exp(-x + a*log(x) - Log_gamma(a)) rather than
-- computing x^a, e^(-x), and Gamma(a) as three separate numbers and
-- combining them after -- the same log-space-first discipline as
-- Log_gamma/Log_beta above, needed for exactly the same reason (x^a and
-- e^(-x) individually over/underflow far before their product does, for
-- large a or x).
local INCOMPLETE_GAMMA_MAX_ITER = 300
local INCOMPLETE_GAMMA_EPS = 1e-15
local INCOMPLETE_GAMMA_TINY = 1e-300

local function gamma_series_p(a, x)
    if x == 0 then
        return 0
    end
    local ap = a
    local total = 1 / a
    local delta = total
    for _ = 1, INCOMPLETE_GAMMA_MAX_ITER do
        ap = ap + 1
        delta = delta * x / ap
        total = total + delta
        if math.abs(delta) < math.abs(total) * INCOMPLETE_GAMMA_EPS then
            break
        end
    end
    return total * math.exp(-x + a * math.log(x) - Log_gamma(a))
end

local function gamma_cf_q(a, x)
    local b = x + 1 - a
    local c = 1 / INCOMPLETE_GAMMA_TINY
    local d = 1 / b
    local h = d
    for i = 1, INCOMPLETE_GAMMA_MAX_ITER do
        local an = -i * (i - a)
        b = b + 2
        d = an * d + b
        if math.abs(d) < INCOMPLETE_GAMMA_TINY then d = INCOMPLETE_GAMMA_TINY end
        c = b + an / c
        if math.abs(c) < INCOMPLETE_GAMMA_TINY then c = INCOMPLETE_GAMMA_TINY end
        d = 1 / d
        local delta = d * c
        h = h * delta
        if math.abs(delta - 1) < INCOMPLETE_GAMMA_EPS then
            break
        end
    end
    return math.exp(-x + a * math.log(x) - Log_gamma(a)) * h
end

-- Gamma_p(a, x) -> P(a,x), the regularized LOWER incomplete gamma
-- function, for a > 0 and x >= 0
--
-- P(a,x) = gamma(a,x)/Gamma(a), where gamma(a,x) = integral from 0 to x
-- of t^(a-1)*e^(-t) dt -- the fraction of Gamma(a)'s defining integral
-- accumulated by the time t reaches x. P(a,0) = 0, P(a,infinity) = 1.
local function Gamma_p(a, x)
    Core.assert_finite_number(a, "Gamma_p: a")
    Core.assert_finite_number(x, "Gamma_p: x")
    assert(a > 0, "Gamma_p: a must be positive.")
    assert(x >= 0, "Gamma_p: x must be nonnegative.")
    if x < a + 1 then
        return gamma_series_p(a, x)
    end
    return 1 - gamma_cf_q(a, x)
end

-- Gamma_q(a, x) -> Q(a,x) = 1 - P(a,x), the regularized UPPER incomplete
-- gamma function, for a > 0 and x >= 0
--
-- Claude: NOT computed as "1 - Gamma_p(a,x)" -- for x large relative to
-- a, P(a,x) itself rounds to exactly 1.0 in double precision (the true
-- Q(a,x) is a genuinely tiny but nonzero number), and "1 - 1.0" would
-- silently give exactly 0 instead. Verified directly: this is exactly
-- what makes Erfc below need this function rather than "1 - Erf(x)" for
-- large x (confirmed in testing against Python: at x=10, "1 - erf(x)"
-- gives exactly 0, while the true erfc(10) is ~2.09e-45).
local function Gamma_q(a, x)
    Core.assert_finite_number(a, "Gamma_q: a")
    Core.assert_finite_number(x, "Gamma_q: x")
    assert(a > 0, "Gamma_q: a must be positive.")
    assert(x >= 0, "Gamma_q: x must be nonnegative.")
    if x < a + 1 then
        return 1 - gamma_series_p(a, x)
    end
    return gamma_cf_q(a, x)
end

-- Erf(x) -> the error function, erf(x) = (2/sqrt(pi)) * integral from 0
-- to x of e^(-t^2) dt
--
-- Built directly on Gamma_p: erf(x) = sign(x) * Gamma_p(1/2, x^2) -- a
-- standard identity (the error function IS a regularized incomplete
-- gamma function in disguise, with a=1/2), reusing the already-verified
-- series/continued-fraction machinery above rather than a separate
-- erf-specific approximation.
local function Erf(x)
    Core.assert_finite_number(x, "Erf: x")
    if x == 0 then
        return 0
    end
    local sign = x > 0 and 1 or -1
    return sign * Gamma_p(0.5, x * x)
end

-- Erfc(x) -> the complementary error function, 1 - Erf(x)
--
-- Claude: NOT computed as "1 - Erf(x)" -- see Gamma_q's own header for
-- why: Erf(x) rounds to exactly 1.0 once x is only moderately large
-- (verified: already happens by x=6), which would silently zero out
-- everything Erfc computes beyond that point, even though the true
-- value keeps shrinking meaningfully for a long time after (verified
-- accurate, against Python's math.erfc, all the way out to x=20, where
-- the true value is ~5.4e-176). Uses Gamma_q(1/2,x^2) directly for
-- x >= 0, and the identity erfc(x) = 2 - erfc(-x) for x < 0.
local function Erfc(x)
    Core.assert_finite_number(x, "Erfc: x")
    if x >= 0 then
        return Gamma_q(0.5, x * x)
    end
    return 2 - Gamma_q(0.5, x * x)
end

-- Bessel_j(n, x) -> J_n(x), the Bessel function of the first kind, order
-- n, for a nonnegative integer n and real x with |x| <= 15
--
-- Claude: computed via the direct power series --
-- J_n(x) = sum_(k=0..infinity) [(-1)^k / (k! * (k+n)!)] * (x/2)^(2k+n)
-- -- which converges for every real x, but whose INTERMEDIATE terms grow
-- before shrinking once x is large enough (roughly x > 2k), causing
-- accumulated rounding error in the alternating sum that gets steadily
-- worse as |x| grows, well before the series stops literally converging.
-- Verified directly against mpmath.besselj (an independent, arbitrary-
-- precision reference) across n = 0..5 and |x| up to 16: relative error
-- stays under 2e-10 for every |x| <= 15, so THAT is the domain enforced
-- here, deliberately conservative -- not the point where the series
-- actually breaks down, but the point up to which this implementation
-- has actually been verified. A large-|x| method (the standard
-- asymptotic expansion, or Miller's backward-recurrence algorithm) is a
-- real, separate future increment -- acknowledged, deferred, not
-- guessed at, same status Bessel_y and the modified Bessel functions
-- have for now.
--
-- J_n(-x) = (-1)^n * J_n(x) falls out of this formula automatically
-- (verified against mpmath), so no separate negative-x branch is needed.
local function Bessel_j(n, x)
    Core.assert_finite_number(x, "Bessel_j: x")
    assert(n == math.floor(n) and n >= 0, "Bessel_j: n must be a nonnegative integer.")
    assert(math.abs(x) <= 15, "Bessel_j: |x| must be <= 15 (the verified domain -- larger |x| is a deferred future increment).")
    if x == 0 then
        return n == 0 and 1 or 0
    end
    local half_x = x / 2
    local factorial_n = 1
    for i = 2, n do factorial_n = factorial_n * i end
    local term = (half_x ^ n) / factorial_n
    local total = term
    local half_x_sq = half_x * half_x
    for k = 1, 200 do
        term = term * (-half_x_sq) / (k * (k + n))
        total = total + term
        if math.abs(term) < math.abs(total) * 1e-17 then
            break
        end
    end
    return total
end

return {
    Gamma = Gamma,
    Log_gamma = Log_gamma,
    Beta = Beta,
    Log_beta = Log_beta,
    Gamma_p = Gamma_p,
    Gamma_q = Gamma_q,
    Erf = Erf,
    Erfc = Erfc,
    Bessel_j = Bessel_j,
}
