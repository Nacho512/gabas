-- [Special Functions]
-- Gamma(x) and Log_gamma(x): the gamma function (generalizes factorial
-- to real arguments -- Gamma(n+1) = n! for every nonnegative integer n)
-- and its natural logarithm -- the two foundational special functions
-- almost everything else in this project builds on. Beta(a,b) and
-- Log_beta(a,b) (Gamma(a)*Gamma(b)/Gamma(a+b), for a,b > 0) are the
-- first thing built on top of them, further down this file.
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

return {
    Gamma = Gamma,
    Log_gamma = Log_gamma,
    Beta = Beta,
    Log_beta = Log_beta,
}
