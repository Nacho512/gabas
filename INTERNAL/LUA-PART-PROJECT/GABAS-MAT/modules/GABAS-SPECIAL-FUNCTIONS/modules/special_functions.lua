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

-- Claude: Beta_i(a, b, x) -> I_x(a,b), the regularized incomplete beta
-- function, for a > 0, b > 0, and x in [0,1]
--
-- I_x(a,b) = B(x;a,b)/B(a,b), where B(x;a,b) = integral from 0 to x of
-- t^(a-1)*(1-t)^(b-1) dt -- the fraction of Beta(a,b)'s defining
-- integral accumulated by the time t reaches x. I_0(a,b) = 0,
-- I_1(a,b) = 1. This is the same continued-fraction technique Gamma_q
-- above already uses (the classical Numerical Recipes "betacf"/"betai"
-- algorithm, transcribed variable-for-variable, not re-derived -- same
-- discipline as the QK21 tables and Brent-Dekker's algorithm), scaled by
-- exp(a*log(x) + b*log(1-x) - Log_beta(a,b)) rather than forming x^a,
-- (1-x)^b, and Beta(a,b) as separate numbers and combining them after --
-- the same log-space-first discipline Log_beta itself already needed,
-- for exactly the same reason.
--
-- The continued fraction converges quickly for x < (a+1)/(a+b+2) and
-- slowly (many more terms, more accumulated rounding error) past it --
-- so past that point, this evaluates it at the SWAPPED, faster-converging
-- point (b, a, 1-x) instead and uses the symmetry identity
-- I_x(a,b) = 1 - I_(1-x)(b,a) to recover the answer, rather than just
-- letting the slow side grind through more iterations. Verified directly
-- against mpmath.betainc (regularized=True) across a wide spread of
-- (a,b,x), including a,b as small as 0.01 and as large as 1000, and x
-- within 1e-6 of either endpoint: relative error stays within ~3.3e-13.
local BETA_I_MAX_ITER = 200
local BETA_I_EPS = 1e-15
local BETA_I_TINY = 1e-300

local function beta_i_cf(a, b, x)
    local qab = a + b
    local qap = a + 1
    local qam = a - 1
    local c = 1.0
    local d = 1 - qab * x / qap
    if math.abs(d) < BETA_I_TINY then d = BETA_I_TINY end
    d = 1 / d
    local h = d
    for m = 1, BETA_I_MAX_ITER do
        local m2 = 2 * m
        local aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1 + aa * d
        if math.abs(d) < BETA_I_TINY then d = BETA_I_TINY end
        c = 1 + aa / c
        if math.abs(c) < BETA_I_TINY then c = BETA_I_TINY end
        d = 1 / d
        h = h * d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1 + aa * d
        if math.abs(d) < BETA_I_TINY then d = BETA_I_TINY end
        c = 1 + aa / c
        if math.abs(c) < BETA_I_TINY then c = BETA_I_TINY end
        d = 1 / d
        local delta = d * c
        h = h * delta
        if math.abs(delta - 1) < BETA_I_EPS then
            break
        end
    end
    return h
end

local function Beta_i(a, b, x)
    Core.assert_finite_number(a, "Beta_i: a")
    Core.assert_finite_number(b, "Beta_i: b")
    Core.assert_finite_number(x, "Beta_i: x")
    assert(a > 0, "Beta_i: a must be positive.")
    assert(b > 0, "Beta_i: b must be positive.")
    assert(x >= 0 and x <= 1, "Beta_i: x must be in [0,1].")
    if x == 0 or x == 1 then
        return x
    end
    local prefactor = math.exp(a * math.log(x) + b * math.log(1 - x) - Log_beta(a, b))
    if x < (a + 1) / (a + b + 2) then
        return prefactor * beta_i_cf(a, b, x) / a
    end
    return 1 - prefactor * beta_i_cf(b, a, 1 - x) / b
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

-- Claude: Bessel_j used to be a direct power series, restricted to
-- |x| <= 15 -- accurate there, but that boundary was a real limitation,
-- not a cosmetic one (the series' intermediate terms grow before
-- shrinking once x is large enough, so accumulated rounding error in the
-- alternating sum gets steadily worse as |x| grows). Replaced here with
-- Miller's algorithm: the classical Numerical Recipes "bessjn"
-- backward-recurrence-plus-normalization technique (independently
-- confirmed, with the same pseudocode, as Algorithm 6/"MILLER_J" in the
-- Bessel reference document Nacho supplied) -- copied/adapted, not
-- derived from scratch, the same discipline as the Lanczos coefficients
-- above.
--
-- The idea: J_n obeys the recurrence J_(k-1)(x) = (2k/x)*J_k(x) -
-- J_(k+1)(x). Run it FORWARD (increasing k) and it is numerically
-- unstable whenever k > x -- rounding error in the direction of the
-- *growing* solution gets amplified every step. Run it BACKWARD
-- (decreasing k, from some arbitrary seed at a high starting order M)
-- and the same recurrence is stable in that direction instead, for any
-- n and x -- the arbitrary seed's error washes out as k decreases. The
-- result is only proportional to the true J_0..J_M, so it still needs
-- normalizing; that uses Neumann's identity
-- 1 = J_0(x) + 2*(J_2(x) + J_4(x) + J_6(x) + ...) -- summing the SAME
-- unnormalized backward-recurrence values, which fixes the scale without
-- needing any separately-computed reference value.
--
-- M must be chosen comfortably above both n and x (too small, and the
-- seed's contamination hasn't washed out by the time the recurrence
-- reaches order n) -- rather than trust a single fixed formula for M,
-- bessel_j_backward_sweep is run at successively larger (always even --
-- required for the Neumann-sum parity to line up correctly, see below)
-- M, starting from max(n, ceil(|x|)) + margin, until the result stops
-- changing -- i.e. self-convergence, no external reference needed at
-- runtime, the same iterate-until-stable pattern gamma_series_p and
-- gamma_cf_q above already use.
--
-- Verified directly against mpmath.besselj (arbitrary-precision
-- reference): relative error stays within ~1e-13 across n = 0..1000 and
-- |x| from 1e-10 to 1000 (including x=0, negative x, and n far larger
-- than x, and x far larger than n) -- a large, systematically swept
-- range, not a handful of lucky points. Genuinely general (no known
-- breakdown mechanism -- unlike the old power series, nothing here
-- degrades as |x| grows, since the recurrence itself never forms a large
-- cancelling sum), but only actually verified up to |x|, n ~ 1000; cost
-- grows with max(n, |x|) (each candidate M costs a full sweep down to 0),
-- so extremely large n or |x| would be slow well before it would be
-- wrong. No upper-bound assertion is enforced -- unlike the old
-- restriction, there is no evidence-based boundary to enforce -- but
-- this hasn't been verified beyond ~1000.
local MILLER_MARGIN = 16
local MILLER_STEP = 16
local MILLER_MAX_TRIES = 100
local MILLER_CONVERGENCE_EPS = 1e-14
local MILLER_RESCALE_BIG = 1e10
local MILLER_RESCALE_SMALL = 1e-10

-- One backward-recurrence sweep from order M down to 0, at fixed
-- positive ax = |x|, returning the corresponding Neumann-normalized
-- J_n(ax). M must be even (see header above).
local function bessel_j_backward_sweep(n, ax, M)
    local tox = 2 / ax
    local bjp = 0 -- J_(j+1), unnormalized
    local bj = 1 -- J_j, unnormalized (arbitrary nonzero seed at j = M)
    local ans = 0
    local sum = 0
    local add_to_sum = false
    for j = M, 1, -1 do
        local bjm = j * tox * bj - bjp
        bjp = bj
        bj = bjm
        if math.abs(bj) > MILLER_RESCALE_BIG then
            bj = bj * MILLER_RESCALE_SMALL
            bjp = bjp * MILLER_RESCALE_SMALL
            ans = ans * MILLER_RESCALE_SMALL
            sum = sum * MILLER_RESCALE_SMALL
        end
        if add_to_sum then
            sum = sum + bj
        end
        add_to_sum = not add_to_sum
        if j == n then
            ans = bjp
        end
    end
    if n == 0 then
        ans = bj -- j never equals n=0 inside the loop (it runs down to 1), so bj after the loop (order 0) is the answer directly.
    end
    sum = 2 * sum - bj
    return ans / sum
end

-- Bessel_j(n, x) -> J_n(x), the Bessel function of the first kind, order
-- n, for a nonnegative integer n and any finite real x
local function Bessel_j(n, x)
    Core.assert_finite_number(x, "Bessel_j: x")
    assert(n == math.floor(n) and n >= 0, "Bessel_j: n must be a nonnegative integer.")
    if x == 0 then
        return n == 0 and 1 or 0
    end
    local ax = math.abs(x)
    local M = math.max(n, math.ceil(ax)) + MILLER_MARGIN
    if M % 2 == 1 then M = M + 1 end
    local value = bessel_j_backward_sweep(n, ax, M)
    for _ = 1, MILLER_MAX_TRIES do
        M = M + MILLER_STEP
        local next_value = bessel_j_backward_sweep(n, ax, M)
        local converged = (next_value == value)
            or math.abs(next_value - value) <= math.abs(next_value) * MILLER_CONVERGENCE_EPS
        value = next_value
        if converged then
            break
        end
    end
    if x < 0 and n % 2 == 1 then
        value = -value
    end
    return value
end

-- Claude: Bessel_y(n,x) -> Y_n(x), the Bessel function of the SECOND
-- kind, order n -- the other independent solution of the same
-- differential equation J_n solves, singular (logarithmically) at
-- x = 0, so only defined here for x > 0.
--
-- Y_0 and Y_1 are computed from their defining series (Abramowitz &
-- Stegun 9.1.11/9.1.13, DLMF 10.8.1) -- copied from that reference, not
-- derived from scratch, same discipline as the Lanczos coefficients
-- above:
--
--   Y_0(x) = (2/pi)*J_0(x)*ln(x/2) - (2/pi) * sum_(k=0..inf) psi(k+1) *
--            (-1)^k * (x/2)^(2k) / (k!)^2
--
--   Y_1(x) = (2/pi)*J_1(x)*ln(x/2) - 2/(pi*x) - (1/pi) *
--            sum_(k=0..inf) (-1)^k * [psi(k+1)+psi(k+2)] *
--            (x/2)^(2k+1) / (k! * (k+1)!)
--
-- where psi(k+1) = -EULER_GAMMA + H_k (H_k the k-th harmonic number,
-- H_0 = 0) is the digamma function at a positive integer -- itself an
-- elementary recurrence, not a separate special function to implement.
--
-- Verified directly against mpmath.bessely: relative error stays under
-- ~5e-12 for 0 < x <= 12, growing past that (same class of hazard as the
-- old Bessel_j power series -- these series' intermediate terms grow
-- before shrinking once x is large enough) -- by x=15 it's already
-- ~5e-11/1e-9, by x=20 ~1e-8. x <= 12 is the domain enforced here,
-- deliberately conservative, exactly the same "verified, not
-- theoretical" discipline the old Bessel_j boundary used. Extending
-- past it (asymptotic expansion for large x, matching how Miller's
-- algorithm extended Bessel_j) is a real, separate future increment.
--
-- Y_n for n >= 2 is built from Y_0, Y_1 by the SAME two-term recurrence
-- Bessel_j uses, run FORWARD (increasing n) this time -- and forward is
-- the numerically stable direction for Y_n (the opposite of J_n): Y_n
-- is the recurrence's dominant, growing solution as n increases, so
-- rounding error in that direction shrinks relative to the true value
-- rather than getting amplified. Verified directly against
-- mpmath.bessely for n up to 20 at x up to 12: no growth in relative
-- error beyond what Y_0/Y_1 already carry at that x.
local EULER_GAMMA = 0.57721566490153286060651209008240243104215933593992
local BESSEL_Y_MAX_X = 12
local BESSEL_Y_SERIES_MAX_ITER = 200
local BESSEL_Y_SERIES_EPS = 1e-18

local function bessel_y0_series(x, j0)
    local half_x_sq = (x / 2) ^ 2
    local term = 1.0
    local harmonic = 0.0
    local total = -EULER_GAMMA * term
    for k = 1, BESSEL_Y_SERIES_MAX_ITER do
        term = term * (-half_x_sq) / (k * k)
        harmonic = harmonic + 1.0 / k
        local contrib = (-EULER_GAMMA + harmonic) * term
        total = total + contrib
        if math.abs(contrib) < BESSEL_Y_SERIES_EPS * math.max(math.abs(total), 1.0) then
            break
        end
    end
    return (2 / math.pi) * j0 * math.log(x / 2) - (2 / math.pi) * total
end

local function bessel_y1_series(x, j1)
    local half_x = x / 2
    local half_x_sq = half_x * half_x
    local fact_k, fact_k1 = 1.0, 1.0
    local harmonic_k, harmonic_k1 = 0.0, 1.0
    local pow_odd = half_x -- (x/2)^(2*0+1)
    local total = (-EULER_GAMMA + (-EULER_GAMMA + harmonic_k1)) * pow_odd / (fact_k * fact_k1)
    local sign = -1
    for k = 1, BESSEL_Y_SERIES_MAX_ITER do
        fact_k = fact_k * k
        fact_k1 = fact_k1 * (k + 1)
        harmonic_k = harmonic_k + 1.0 / k
        harmonic_k1 = harmonic_k1 + 1.0 / (k + 1)
        pow_odd = pow_odd * half_x_sq
        local psi_sum = (-EULER_GAMMA + harmonic_k) + (-EULER_GAMMA + harmonic_k1)
        local contrib = sign * psi_sum * pow_odd / (fact_k * fact_k1)
        total = total + contrib
        sign = -sign
        if math.abs(contrib) < BESSEL_Y_SERIES_EPS * math.max(math.abs(total), 1.0) then
            break
        end
    end
    return (2 / math.pi) * j1 * math.log(x / 2) - 2 / (math.pi * x) - (1 / math.pi) * total
end

local function Bessel_y(n, x)
    Core.assert_finite_number(x, "Bessel_y: x")
    assert(n == math.floor(n) and n >= 0, "Bessel_y: n must be a nonnegative integer.")
    assert(x > 0, "Bessel_y: x must be positive (Y_n has a logarithmic singularity at x=0 and is not real-valued for x<0).")
    assert(x <= BESSEL_Y_MAX_X, "Bessel_y: x must be <= 12 (the verified domain -- larger x is a deferred future increment).")
    local y0 = bessel_y0_series(x, Bessel_j(0, x))
    if n == 0 then
        return y0
    end
    local y1 = bessel_y1_series(x, Bessel_j(1, x))
    if n == 1 then
        return y1
    end
    local y_prev, y_curr = y0, y1
    for k = 1, n - 1 do
        local y_next = (2 * k / x) * y_curr - y_prev
        y_prev, y_curr = y_curr, y_next
    end
    return y_curr
end

-- Modified_bessel_i(n, x) -> I_n(x), the modified Bessel function of the
-- first kind, order n, for a nonnegative integer n and any finite real x
--
-- Claude: I_n(x) = sum_(k=0..infinity) [(x/2)^(2k+n)] / (k! * (k+n)!) --
-- the same defining structure as Bessel_j's power series, but WITHOUT
-- the alternating (-1)^k sign. Every term is the same sign (nonnegative
-- for x >= 0), so unlike Bessel_j's series there is no cancellation at
-- all -- nothing to accumulate rounding error, at any x. Verified
-- directly against mpmath.besseli across n = 0..50 and x up to 750:
-- relative error stays at the few-1e-16 (machine epsilon) level
-- everywhere the true result is representable as a double, with no
-- degradation as x grows (unlike Bessel_j's old power series) -- the
-- only failure mode is I_n(x) itself exceeding a double's range, which
-- happens above roughly x=709 (depending on n) and returns +math.huge,
-- the same honest IEEE overflow Gamma(x) already has beyond x~172, not
-- a silently wrong finite number. No domain restriction is enforced
-- here for exactly that reason -- there is no accuracy boundary to
-- enforce, only the same overflow every double-precision routine has.
--
-- A separate EXPONENTIALLY-SCALED variant (I_n(x)*e^-x, representable
-- even where I_n(x) itself overflows -- see Log_gamma's own header for
-- the same idea applied to Gamma) is deferred, not needed by anything
-- in this file yet -- but the Bessel reference document Nacho supplied
-- flags exactly this as a required building block for
-- Modified_bessel_k, so it will come back when that's built.
--
-- I_n(-x) = (-1)^n * I_n(x) falls out of this formula automatically
-- (the sign of half_x^n for odd n is handled explicitly below, since
-- the log-space computation that produces it needs abs(half_x) --
-- log() of a negative number isn't real), so no separate negative-x
-- CASE needs its own algorithm -- same as Bessel_j.
--
-- Claude: the very first series term, (half_x^n)/n!, is computed in
-- LOG-SPACE FIRST (n*log(|half_x|) - Log_gamma(n+1), then exp() of
-- that) -- not as half_x^n and n! formed as separate numbers and then
-- divided -- the same log-space-first discipline as Gamma/Log_gamma
-- above, needed for exactly the same reason: half_x^n alone can
-- overflow a double even when the true ratio (half_x^n)/n! is
-- comfortably representable. Caught directly by testing, not assumed:
-- Modified_bessel_i(150, 300) computed the naive way returns
-- +math.huge (half_x^150 = 150^150 overflows on its own), while the
-- true value is an ordinary, representable ~4.5e112 -- confirmed
-- against mpmath before this fix. An earlier version of this function
-- also computed n! via a Lua-integer accumulator (`factorial_n = 1`
-- rather than `1.0`) -- Lua 5.4 integers are 64-bit and WRAP SILENTLY
-- on overflow instead of promoting to float, and 50! already exceeds
-- that range -- also caught directly by testing
-- (Modified_bessel_i(50, 100) came back negative). Routing n! through
-- Log_gamma(n+1) (already float-based throughout) sidesteps that
-- failure mode entirely, not just patches it.
local MODIFIED_BESSEL_I_MAX_ITER = 600

local function Modified_bessel_i(n, x)
    Core.assert_finite_number(x, "Modified_bessel_i: x")
    assert(n == math.floor(n) and n >= 0, "Modified_bessel_i: n must be a nonnegative integer.")
    if x == 0 then
        return n == 0 and 1 or 0
    end
    local half_x = x / 2
    local abs_half_x = math.abs(half_x)
    local sign0 = (half_x < 0 and n % 2 == 1) and -1 or 1
    local term = sign0 * math.exp(n * math.log(abs_half_x) - Log_gamma(n + 1))
    local total = term
    local half_x_sq = half_x * half_x
    for k = 1, MODIFIED_BESSEL_I_MAX_ITER do
        term = term * half_x_sq / (k * (k + n))
        total = total + term
        if total == math.huge then
            break
        end
        if math.abs(term) < math.abs(total) * 1e-17 then
            break
        end
    end
    return total
end

-- Claude: Modified_bessel_k(n,x) -> K_n(x), the modified Bessel function
-- of the SECOND kind, order n -- the other independent solution of the
-- same ODE Modified_bessel_i solves, singular at x = 0, so only defined
-- here for x > 0. This is the hardest of the four Bessel families the
-- reference document Nacho supplied covers (Section 19 flags dedicated
-- near-integer-order methods as needing a verified publication, not
-- casual reconstruction) -- but that specific hazard is about
-- NONINTEGER order, where the connection formula
-- K_nu = (pi/2)*(I_(-nu)-I_nu)/sin(pi*nu) becomes ill-conditioned near
-- integer nu. Restricted to nonnegative INTEGER n like every other
-- Bessel function in this file, that whole problem doesn't arise --
-- K_0 and K_1 are computed directly from their own defining series
-- (Abramowitz & Stegun 9.6.13/9.6.15), no connection formula involved.
--
--   K_0(x) = -[ln(x/2)+EULER_GAMMA]*I_0(x) + sum_(k=1..inf) H_k *
--            (x/2)^(2k) / (k!)^2
--
--   K_1(x) = 1/x + [ln(x/2)+EULER_GAMMA]*I_1(x) - (x/4) *
--            sum_(k=0..inf) (H_k+H_(k+1)) / (k!*(k+1)!) * (x/2)^(2k)
--
-- (H_k the k-th harmonic number, same building block Bessel_y already
-- uses -- psi(k+1) = -EULER_GAMMA + H_k, so these are the same kind of
-- series, just without the alternating sign Y's has, since I_n's series
-- (unlike J_n's) has none either).
--
-- Verified directly against mpmath.besselk: relative error stays under
-- ~7e-12 for 0 < x <= 6 -- a MUCH narrower domain than Bessel_y's x<=12
-- despite the strong structural similarity of the two series, because
-- I_0(x)/I_1(x) here GROW exponentially while the true K_0(x)/K_1(x)
-- DECAY exponentially: the (ln(x/2)+EULER_GAMMA)*I_n(x) term and the
-- series total end up catastrophically cancelling against each other,
-- and that cancellation gets worse far faster than Bessel_y's analogous
-- (much milder) cancellation against a bounded J_n(x). By x=10 the
-- error is already ~1e-7. x <= 6 is deliberately conservative -- same
-- "verified, not theoretical" discipline as every other domain boundary
-- in this file. Extending past it needs the large-x asymptotic
-- expansion (in terms of exponentially SCALED I_nu/K_nu, per the
-- reference document -- unscaled K_nu underflows long before an
-- asymptotic expansion would even be needed) -- deferred, real future
-- work, not attempted here.
--
-- K_n for n >= 2 is built from K_0, K_1 by the same forward recurrence
-- style as Bessel_y (forward is stable for K_n too, for the same
-- reason: K_n is the dominant, growing solution as n increases) --
-- verified directly against mpmath.besselk for n up to 20 at x up to 6:
-- no growth in relative error beyond what K_0/K_1 already carry. Note
-- the PLUS sign here (K_(n+1) = K_(n-1) + (2n/x)*K_n), unlike
-- Bessel_j/Bessel_y's minus -- K's defining ODE has the opposite sign
-- on its x^2 term, which flips this recurrence's sign too.
local MODIFIED_BESSEL_K_MAX_X = 6
local MODIFIED_BESSEL_K_SERIES_MAX_ITER = 300
local MODIFIED_BESSEL_K_SERIES_EPS = 1e-18

local function bessel_k0_series(x, i0)
    local half_x_sq = (x / 2) ^ 2
    local term = 1.0
    local harmonic = 0.0
    local total = 0.0
    for k = 1, MODIFIED_BESSEL_K_SERIES_MAX_ITER do
        term = term * half_x_sq / (k * k)
        harmonic = harmonic + 1.0 / k
        local contrib = harmonic * term
        total = total + contrib
        if math.abs(contrib) < MODIFIED_BESSEL_K_SERIES_EPS * math.max(math.abs(total), 1.0) then
            break
        end
    end
    return -(math.log(x / 2) + EULER_GAMMA) * i0 + total
end

local function bessel_k1_series(x, i1)
    local half_x = x / 2
    local half_x_sq = half_x * half_x
    local fact_k, fact_k1 = 1.0, 1.0
    local harmonic_k, harmonic_k1 = 0.0, 1.0
    local pow_even = 1.0 -- (x/2)^(2*0)
    local total = (harmonic_k + harmonic_k1) * pow_even / (fact_k * fact_k1) -- k=0 term
    for k = 1, MODIFIED_BESSEL_K_SERIES_MAX_ITER do
        fact_k = fact_k * k
        fact_k1 = fact_k1 * (k + 1)
        harmonic_k = harmonic_k + 1.0 / k
        harmonic_k1 = harmonic_k1 + 1.0 / (k + 1)
        pow_even = pow_even * half_x_sq
        local contrib = (harmonic_k + harmonic_k1) * pow_even / (fact_k * fact_k1)
        total = total + contrib
        if math.abs(contrib) < MODIFIED_BESSEL_K_SERIES_EPS * math.max(math.abs(total), 1.0) then
            break
        end
    end
    return 1 / x + (math.log(x / 2) + EULER_GAMMA) * i1 - (half_x / 2) * total
end

local function Modified_bessel_k(n, x)
    Core.assert_finite_number(x, "Modified_bessel_k: x")
    assert(n == math.floor(n) and n >= 0, "Modified_bessel_k: n must be a nonnegative integer.")
    assert(x > 0, "Modified_bessel_k: x must be positive (K_n has a singularity at x=0 and is not real-valued for x<0).")
    assert(x <= MODIFIED_BESSEL_K_MAX_X, "Modified_bessel_k: x must be <= 6 (the verified domain -- larger x is a deferred future increment).")
    local k0 = bessel_k0_series(x, Modified_bessel_i(0, x))
    if n == 0 then
        return k0
    end
    local k1 = bessel_k1_series(x, Modified_bessel_i(1, x))
    if n == 1 then
        return k1
    end
    local k_prev, k_curr = k0, k1
    for j = 1, n - 1 do
        local k_next = k_prev + (2 * j / x) * k_curr
        k_prev, k_curr = k_curr, k_next
    end
    return k_curr
end

return {
    Gamma = Gamma,
    Log_gamma = Log_gamma,
    Beta = Beta,
    Log_beta = Log_beta,
    Beta_i = Beta_i,
    Gamma_p = Gamma_p,
    Gamma_q = Gamma_q,
    Erf = Erf,
    Erfc = Erfc,
    Bessel_j = Bessel_j,
    Bessel_y = Bessel_y,
    Modified_bessel_i = Modified_bessel_i,
    Modified_bessel_k = Modified_bessel_k,
}
