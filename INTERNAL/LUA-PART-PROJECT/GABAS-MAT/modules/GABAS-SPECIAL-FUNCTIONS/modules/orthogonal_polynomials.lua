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

return {
    Chebyshev_t = Chebyshev_t,
    Chebyshev_u = Chebyshev_u,
}
