-- [Orthogonal Polynomials]
-- Chebyshev_t(n,x) -> T_n(x) and Chebyshev_u(n,x) -> U_n(x), the
-- Chebyshev polynomials of the first and second kind. Split into its
-- own file (rather than folded into special_functions.lua) since this
-- is a genuinely different family -- polynomials built from simple
-- three-term recurrences, not the series/continued-fraction/asymptotic
-- machinery the rest of GABAS-SPECIAL-FUNCTIONS needs -- and more
-- orthogonal-polynomial families (Legendre, Hermite, Laguerre, ...) are
-- planned to join this file incrementally.
local Core = require("GABAS-SPECIAL-FUNCTIONS.modules.core")

-- Chebyshev_t(n, x) -> T_n(x), the Chebyshev polynomial of the first
-- kind, order n, for a nonnegative integer n and any finite real x.
--
-- Claude: computed via the direct three-term recurrence
-- T_(k+1)(x) = 2x*T_k(x) - T_(k-1)(x), T_0=1, T_1=x -- unlike
-- Bessel_j's own three-term recurrence, this one needs no Miller's-
-- algorithm trick: forward recurrence here computes the DOMINANT
-- solution in every regime -- for |x|>1, T_n genuinely grows like x^n
-- (that's the correct, expected behavior, not overflow-driven error);
-- for |x|<=1, T_n stays bounded (|T_n(x)|<=1) and forward iteration is
-- neutrally stable (a rounding error introduced at one step doesn't get
-- amplified by later ones, only carried forward at the same
-- magnitude). This is the standard, well-known technique -- no external
-- reference needed beyond directly testing it, same as this project's
-- other simple textbook algorithms (factorial, GCD, bisection, ...).
--
-- Verified directly against the independent closed forms
-- T_n(x) = cos(n*arccos(x)) for |x|<=1 and T_n(x) = cosh(n*arccosh(x))
-- for x>1 (T_n(-x) = (-1)^n*T_n(x) covers x<-1) -- relative error stays
-- at the few-1e-15 (machine epsilon) level for n up to 100 and |x| up
-- to 5.
--
-- Claude: "x = x + 0.0" up front is deliberate, not a style nit -- the
-- same Lua-integer-overflow hazard Modified_bessel_i's own header
-- documents. A caller passing an integer-looking argument (e.g.
-- Chebyshev_t(34, -2), no decimal point) hands x in as a Lua *integer*,
-- and 1 (the seed) is one too -- so without this line, the entire
-- recurrence would run in 64-bit integer arithmetic and wrap silently
-- once |T_n(x)| exceeds ~9.2e18 (T_34(-2) already does), rather than
-- promoting to a float that can represent it. Caught directly by
-- testing Chebyshev_t(34, -2) against an independent reference (it
-- came back negative and wildly wrong), not assumed.
local function Chebyshev_t(n, x)
    Core.assert_finite_number(x, "Chebyshev_t: x")
    assert(n == math.floor(n) and n >= 0, "Chebyshev_t: n must be a nonnegative integer.")
    x = x + 0.0
    if n == 0 then
        return 1.0
    end
    local t_prev, t_curr = 1.0, x
    for _ = 2, n do
        t_prev, t_curr = t_curr, 2 * x * t_curr - t_prev
    end
    return t_curr
end

-- Chebyshev_u(n, x) -> U_n(x), the Chebyshev polynomial of the second
-- kind, order n, for a nonnegative integer n and any finite real x.
--
-- Claude: the same recurrence family as Chebyshev_t above, just a
-- different seed (U_0=1, U_1=2x, then the identical
-- U_(k+1)(x) = 2x*U_k(x) - U_(k-1)(x)) -- same stability reasoning
-- applies here too.
--
-- Verified directly against the independent closed forms
-- U_n(x) = sin((n+1)*arccos(x))/sin(arccos(x)) for |x|<1 and the
-- analogous sinh/arccosh form for |x|>1 -- relative error stays under
-- ~5e-14 for n up to 100 and |x| up to 5.
local function Chebyshev_u(n, x)
    Core.assert_finite_number(x, "Chebyshev_u: x")
    assert(n == math.floor(n) and n >= 0, "Chebyshev_u: n must be a nonnegative integer.")
    x = x + 0.0 -- Claude: force float -- same reasoning as Chebyshev_t above.
    if n == 0 then
        return 1.0
    end
    local u_prev, u_curr = 1.0, 2 * x
    for _ = 2, n do
        u_prev, u_curr = u_curr, 2 * x * u_curr - u_prev
    end
    return u_curr
end

-- Legendre_p(n, x) -> P_n(x), the Legendre polynomial of order n, for a
-- nonnegative integer n and any finite real x.
--
-- Claude: computed via the direct three-term recurrence
-- (k+1)*P_(k+1)(x) = (2k+1)*x*P_k(x) - k*P_(k-1)(x), P_0=1, P_1=x --
-- the standard, well-known technique (Bonnet's recursion formula), same
-- "no external reference needed, just test it directly" status as
-- Chebyshev_t/Chebyshev_u above. Unlike those, this recurrence divides
-- by (k+1) at each step rather than only ever doubling and subtracting
-- -- still exact, well-conditioned arithmetic (no cancellation: for
-- |x|<=1, the true P_n(x) stays bounded in [-1,1] and every
-- intermediate value along the way does too; for |x|>1, P_n(x) grows,
-- again the genuine dominant solution, not overflow-driven error).
--
-- Verified directly against mpmath.legendre (an independent, trusted
-- reference): relative error stays within ~6e-15 for n up to 100 and
-- |x| up to 5.
local function Legendre_p(n, x)
    Core.assert_finite_number(x, "Legendre_p: x")
    assert(n == math.floor(n) and n >= 0, "Legendre_p: n must be a nonnegative integer.")
    -- Claude: "x + 0.0" -- same defensive habit as Chebyshev_t/
    -- Chebyshev_u, for consistency -- but NOT load-bearing here the way
    -- it was there: Lua's "/" operator always returns a float even for
    -- two integer operands, and this recurrence divides by (k+1) every
    -- single step, so it self-promotes to float arithmetic on its very
    -- first iteration regardless. Checked directly, not assumed:
    -- Legendre_p(34, -2) (an integer-only call, deliberately mirroring
    -- the Chebyshev_t regression case) already comes back correct
    -- (~2.797e18, matching mpmath.legendre) even without this line.
    x = x + 0.0
    if n == 0 then
        return 1.0
    end
    local p_prev, p_curr = 1.0, x
    for k = 1, n - 1 do
        p_prev, p_curr = p_curr, ((2 * k + 1) * x * p_curr - k * p_prev) / (k + 1)
    end
    return p_curr
end

-- Hermite_h(n, x) -> H_n(x), the (physicists') Hermite polynomial of
-- order n, for a nonnegative integer n and any finite real x.
--
-- Claude: computed via the direct three-term recurrence
-- H_(k+1)(x) = 2x*H_k(x) - 2k*H_(k-1)(x), H_0=1, H_1=2x -- same
-- "well-known, test it directly" status as Chebyshev_t/Chebyshev_u/
-- Legendre_p above. H_n(x) grows fast (like (2x)^n, and the 2k
-- coefficient grows too) -- genuine growth, not overflow-driven error,
-- and Gamma already overflows to +-math.huge past a similar order of
-- magnitude in this same file, so nothing new is needed here for that.
--
-- Verified directly against mpmath.hermite (an independent, trusted
-- reference): relative error stays within ~9e-16 for n up to 50 and
-- |x| up to 5.
local function Hermite_h(n, x)
    Core.assert_finite_number(x, "Hermite_h: x")
    assert(n == math.floor(n) and n >= 0, "Hermite_h: n must be a nonnegative integer.")
    x = x + 0.0 -- Claude: force float -- same defensive habit as Chebyshev_t/Chebyshev_u/Legendre_p above.
    if n == 0 then
        return 1.0
    end
    local h_prev, h_curr = 1.0, 2 * x
    for k = 1, n - 1 do
        h_prev, h_curr = h_curr, 2 * x * h_curr - 2 * k * h_prev
    end
    return h_curr
end

-- Laguerre_l(n, x) -> L_n(x), the (plain, alpha=0) Laguerre polynomial
-- of order n, for a nonnegative integer n and any finite real x.
--
-- Claude: computed via the direct three-term recurrence
-- (k+1)*L_(k+1)(x) = (2k+1-x)*L_k(x) - k*L_(k-1)(x), L_0=1, L_1=1-x --
-- same "well-known, test it directly" status as the other orthogonal
-- polynomials in this file. Like Legendre_p, this divides by (k+1)
-- every step, so it self-promotes to float arithmetic on its first
-- iteration regardless of whether x was passed as a Lua integer --
-- checked directly (Laguerre_l(34, -2) stays accurate without needing
-- explicit float coercion), same situation as Legendre_p, for the same
-- reason. Kept the "x + 0.0" line anyway for consistency across this
-- file.
--
-- Verified directly against mpmath.laguerre(n, 0, x) (an independent,
-- trusted reference): relative error stays within ~4.5e-14 for n up to
-- 80 and x from -2 to 10 -- unlike the other polynomials in this file,
-- Laguerre's natural domain of interest is x >= 0 (it's the family
-- orthogonal on [0, infinity) with weight e^-x), but the recurrence
-- itself is a perfectly well-defined polynomial identity for any real
-- x, so no domain restriction is imposed here, matching how
-- Chebyshev_t/Chebyshev_u/Legendre_p/Hermite_h don't restrict to their
-- own "natural" intervals either.
--
-- The generalized (associated) form L_n^alpha(x), needed for e.g. the
-- radial part of the hydrogen atom's wavefunctions, is a real, separate
-- future increment -- not attempted here (alpha=0 is the special case
-- this implements).
local function Laguerre_l(n, x)
    Core.assert_finite_number(x, "Laguerre_l: x")
    assert(n == math.floor(n) and n >= 0, "Laguerre_l: n must be a nonnegative integer.")
    x = x + 0.0
    if n == 0 then
        return 1.0
    end
    local l_prev, l_curr = 1.0, 1 - x
    for k = 1, n - 1 do
        l_prev, l_curr = l_curr, ((2 * k + 1 - x) * l_curr - k * l_prev) / (k + 1)
    end
    return l_curr
end

return {
    Chebyshev_t = Chebyshev_t,
    Chebyshev_u = Chebyshev_u,
    Legendre_p = Legendre_p,
    Hermite_h = Hermite_h,
    Laguerre_l = Laguerre_l,
}
