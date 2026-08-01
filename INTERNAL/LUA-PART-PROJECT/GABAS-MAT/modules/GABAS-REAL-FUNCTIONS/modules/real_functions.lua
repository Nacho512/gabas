-- [Real Functions]
-- The trivial, genuinely-new real one-variable utility functions this
-- project needs: Sign, Floor, Ceil, Round, Trunc, Frac, Min, Max, Clamp,
-- Hypot.
--
-- Claude: deliberately does NOT include sin/cos/tan/asin/acos/atan/
-- sinh/cosh/tanh/exp/ln/sqrt/abs -- those already exist, properly, as
-- Dual/Complex-aware dispatchers inside
-- GABAS-CALC-ONE-VAR/modules/compile_expression.lua (needed there for
-- Derivative_at_point to differentiate through them). Duplicating them
-- here as real-only wrappers would create two "sin" functions with
-- different capabilities living in two projects -- exactly the
-- duplication risk this whole codebase has been careful to avoid.
-- ARCHITECTURALLY, those thirteen belong in this project too, and
-- compile_expression.lua should eventually depend on this module rather
-- than define its own copies -- but that migration is its own careful
-- future step, not folded into this one.
--
-- Every function here is REAL-only (rejects Complex and Dual values
-- outright, via Core.assert_finite_number) -- not an arbitrary
-- restriction: Sign/Floor/Ceil/Round/Trunc/Frac/Min/Max/Clamp are all
-- either piecewise-constant or have jump discontinuities/kinks, so
-- "differentiate through this with Dual numbers" would not even be
-- mathematically meaningful for most of them. Hypot alone is smooth
-- away from the origin, but is kept real-only here too for consistency
-- -- a Dual-aware Hypot is a real future increment if/when something
-- actually needs to differentiate through it, not added speculatively.
local Core = require("GABAS-REAL-FUNCTIONS.modules.core")

-- Sign(x) -> -1, 0, or 1 -- 0 by convention at x = 0 (not "no sign" or
-- an error; the standard mathematical convention that keeps sign(0)*0 =
-- 0 consistent with sign(x)*|x| = x everywhere).
local function Sign(x)
    Core.assert_finite_number(x, "Sign: x")
    if x > 0 then return 1 end
    if x < 0 then return -1 end
    return 0
end

-- Floor(x) -> the largest integer <= x. A thin, validated wrapper around
-- Lua's own math.floor -- Lua already provides this natively and
-- correctly; the value added here is this project's own input contract
-- (a clear, specific error for a non-finite x), not a new algorithm.
local function Floor(x)
    Core.assert_finite_number(x, "Floor: x")
    return math.floor(x)
end

-- Ceil(x) -> the smallest integer >= x. Same reasoning as Floor.
local function Ceil(x)
    Core.assert_finite_number(x, "Ceil: x")
    return math.ceil(x)
end

-- Round(x) -> the nearest integer to x, ties rounding AWAY from zero
-- (2.5 -> 3, -2.5 -> -3) -- one of several defensible tie-breaking
-- conventions (round-half-to-even is another common one); this is the
-- one most calculators and everyday arithmetic use, and the one chosen
-- here, stated explicitly so it is never ambiguous which convention a
-- caller is getting.
local function Round(x)
    Core.assert_finite_number(x, "Round: x")
    if x >= 0 then
        return math.floor(x + 0.5)
    end
    return math.ceil(x - 0.5)
end

-- Trunc(x) -> x with any fractional part discarded, rounding TOWARD
-- zero (Trunc(2.7) = 2, Trunc(-2.7) = -2) -- distinct from Floor, which
-- always rounds DOWN regardless of sign (Floor(-2.7) = -3).
local function Trunc(x)
    Core.assert_finite_number(x, "Trunc: x")
    if x >= 0 then
        return math.floor(x)
    end
    return math.ceil(x)
end

-- Frac(x) -> the fractional part of x, x - Trunc(x) -- always the same
-- sign as x (or zero), by construction (Frac(2.7) = 0.7,
-- Frac(-2.7) = -0.7).
local function Frac(x)
    Core.assert_finite_number(x, "Frac: x")
    return x - Trunc(x)
end

-- Min(a, b) -> the smaller of a and b.
local function Min(a, b)
    Core.assert_finite_number(a, "Min: a")
    Core.assert_finite_number(b, "Min: b")
    if a < b then return a end
    return b
end

-- Max(a, b) -> the larger of a and b.
local function Max(a, b)
    Core.assert_finite_number(a, "Max: a")
    Core.assert_finite_number(b, "Max: b")
    if a > b then return a end
    return b
end

-- Clamp(x, lower, upper) -> x restricted to the closed interval
-- [lower, upper] -- x itself if it already lies inside, otherwise
-- whichever bound it crossed.
local function Clamp(x, lower, upper)
    Core.assert_finite_number(x, "Clamp: x")
    Core.assert_finite_number(lower, "Clamp: lower")
    Core.assert_finite_number(upper, "Clamp: upper")
    assert(lower <= upper, "Clamp: lower must not exceed upper.")
    if x < lower then return lower end
    if x > upper then return upper end
    return x
end

-- Hypot(x, y) -> sqrt(x^2 + y^2), computed via the standard OVERFLOW-
-- SAFE algorithm rather than literally squaring and adding -- naive
-- sqrt(x*x + y*y) overflows to +inf whenever |x| or |y| exceeds roughly
-- sqrt(the largest representable double) (~1.3e154), even though the
-- true hypotenuse is still a perfectly ordinary, representable number
-- (verified directly: Hypot(1e200, 1e200) below returns a finite value,
-- ~1.414e200, while the naive formula would silently return +inf).
--
-- Claude: factors out the larger magnitude first -- ax*sqrt(1+(ay/ax)^2)
-- with ax = max(|x|,|y|) -- so the squared term (ay/ax)^2 is always in
-- [0,1] and can never overflow, and the final multiplication by ax never
-- introduces overflow beyond what ax itself, one of the original inputs,
-- already represented safely. This is the same core technique (rescale
-- before squaring, don't square first) as the safe_hypot/overflow-safe
-- norm computation referenced throughout this project's own numerical-
-- rigor discipline -- well-known, not improvised.
local function Hypot(x, y)
    Core.assert_finite_number(x, "Hypot: x")
    Core.assert_finite_number(y, "Hypot: y")
    local ax, ay = math.abs(x), math.abs(y)
    if ax < ay then ax, ay = ay, ax end
    if ax == 0 then return 0 end
    local ratio = ay / ax
    return ax * math.sqrt(1 + ratio * ratio)
end

return {
    Sign = Sign,
    Floor = Floor,
    Ceil = Ceil,
    Round = Round,
    Trunc = Trunc,
    Frac = Frac,
    Min = Min,
    Max = Max,
    Clamp = Clamp,
    Hypot = Hypot,
}
