local RealFunctions = require("GABAS-REAL-FUNCTIONS.GABAS_REAL_FUNCTIONS")
local Testing = require("GABAS-LINAL.modules.testing")

-- ===== Sign =====

assert(RealFunctions.Sign(5) == 1, "Sign(5)")
assert(RealFunctions.Sign(-5) == -1, "Sign(-5)")
assert(RealFunctions.Sign(0) == 0, "Sign(0) == 0 by convention")
Testing.Assert_error(function() return RealFunctions.Sign("5") end, "Sign: x")

-- ===== Floor / Ceil =====

assert(RealFunctions.Floor(2.7) == 2, "Floor(2.7)")
assert(RealFunctions.Floor(-2.7) == -3, "Floor(-2.7) rounds DOWN regardless of sign")
assert(RealFunctions.Ceil(2.1) == 3, "Ceil(2.1)")
assert(RealFunctions.Ceil(-2.1) == -2, "Ceil(-2.1)")

-- ===== Round (ties away from zero) =====

assert(RealFunctions.Round(2.5) == 3, "Round(2.5), tie away from zero")
assert(RealFunctions.Round(-2.5) == -3, "Round(-2.5), tie away from zero")
assert(RealFunctions.Round(2.4) == 2, "Round(2.4)")
assert(RealFunctions.Round(2.6) == 3, "Round(2.6)")

-- ===== Trunc (distinct from Floor for negative x) =====

assert(RealFunctions.Trunc(2.7) == 2, "Trunc(2.7)")
assert(RealFunctions.Trunc(-2.7) == -2, "Trunc(-2.7) rounds TOWARD zero, unlike Floor")

-- ===== Frac =====

Testing.Assert_close(RealFunctions.Frac(2.7), 0.7, 1e-9, "Frac(2.7)")
Testing.Assert_close(RealFunctions.Frac(-2.7), -0.7, 1e-9, "Frac(-2.7), same sign as x")
assert(RealFunctions.Frac(5) == 0, "Frac of an integer is 0")

-- ===== Min / Max =====

assert(RealFunctions.Min(3, 7) == 3, "Min(3,7)")
assert(RealFunctions.Min(7, 3) == 3, "Min(7,3), order-independent")
assert(RealFunctions.Max(3, 7) == 7, "Max(3,7)")
assert(RealFunctions.Max(-3, -7) == -3, "Max(-3,-7)")

-- ===== Clamp =====

assert(RealFunctions.Clamp(5, 0, 10) == 5, "Clamp(5,0,10), already inside")
assert(RealFunctions.Clamp(-5, 0, 10) == 0, "Clamp(-5,0,10), below lower bound")
assert(RealFunctions.Clamp(15, 0, 10) == 10, "Clamp(15,0,10), above upper bound")
Testing.Assert_error(function() return RealFunctions.Clamp(5, 10, 0) end, "Clamp: lower must not exceed upper")

-- ===== Hypot =====

Testing.Assert_close(RealFunctions.Hypot(3, 4), 5, 1e-12, "Hypot(3,4), the 3-4-5 triangle")
assert(RealFunctions.Hypot(0, 0) == 0, "Hypot(0,0)")
Testing.Assert_close(RealFunctions.Hypot(-3, 4), 5, 1e-12, "Hypot is sign-independent")

-- The overflow-safe algorithm: naive sqrt(x*x+y*y) would compute
-- x*x = 1e400, which overflows a double (max ~1.8e308) to +inf, giving
-- a WRONG answer of +inf for a perfectly ordinary, representable true
-- hypotenuse (~1.41421e200). Verified against an independent
-- high-precision reference (Python's math.hypot, which uses the same
-- overflow-safe technique).
local big = RealFunctions.Hypot(1e200, 1e200)
assert(big ~= math.huge, "Hypot must not overflow for large but representable inputs")
Testing.Assert_close(big, 1.4142135623730951e200, 1e190, "Hypot(1e200,1e200) matches the true value")

-- Confirm the NAIVE formula really would overflow here, so the test
-- above is actually exercising the overflow-safety, not a redundant
-- check.
assert(1e200 * 1e200 == math.huge, "sanity check: naive x*x really does overflow at this scale")

Testing.Assert_error(function() return RealFunctions.Hypot("3", 4) end, "Hypot: x")
