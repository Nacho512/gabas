local SpecialFunctions = require("GABAS-SPECIAL-FUNCTIONS.GABAS_SPECIAL_FUNCTIONS")
local Combinatorics = require("GABAS-COMBINATORICS.GABAS_COMBINATORICS")
local Testing = require("GABAS-LINAL.modules.testing")

-- ===== Gamma: known closed-form and integer values =====

-- Gamma(n+1) = n! for every nonnegative integer n -- cross-checked
-- directly against GABAS-COMBINATORICS's own, independently verified
-- Factorial (differential testing, not just internal self-consistency).
for n = 0, 10 do
    Testing.Assert_close(SpecialFunctions.Gamma(n + 1), Combinatorics.Factorial(n), 1e-6,
        "Gamma(" .. (n + 1) .. ") == " .. n .. "!")
end

-- Gamma(1/2) = sqrt(pi) (the classical closed form).
Testing.Assert_close(SpecialFunctions.Gamma(0.5), math.sqrt(math.pi), 1e-12, "Gamma(0.5) = sqrt(pi)")

-- Gamma(-1/2) = -2*sqrt(pi) -- exercises the reflection formula AND its
-- sign handling for negative x.
Testing.Assert_close(SpecialFunctions.Gamma(-0.5), -2 * math.sqrt(math.pi), 1e-12, "Gamma(-0.5) = -2*sqrt(pi)")

-- Reference values below (mpmath-independent: Python's own math.gamma)
-- for non-integer, non-half-integer x.
Testing.Assert_close(SpecialFunctions.Gamma(0.3), 2.991568987687591, 1e-9, "Gamma(0.3)")
Testing.Assert_close(SpecialFunctions.Gamma(3.7), 4.170651783796603, 1e-9, "Gamma(3.7)")
Testing.Assert_close(SpecialFunctions.Gamma(-4.5), -0.06001960130050425, 1e-9, "Gamma(-4.5)")

-- Poles at every nonpositive integer -- must be rejected explicitly, not
-- silently returned as Inf/NaN.
Testing.Assert_error(function() return SpecialFunctions.Gamma(0) end, "Gamma: x must not be a nonpositive integer")
Testing.Assert_error(function() return SpecialFunctions.Gamma(-1) end, "Gamma: x must not be a nonpositive integer")
Testing.Assert_error(function() return SpecialFunctions.Gamma(-7) end, "Gamma: x must not be a nonpositive integer")
Testing.Assert_error(function() return SpecialFunctions.Gamma("2") end, "Gamma: x")

-- ===== Log_gamma =====
--
-- Reference values below are Python's math.lgamma (an independent,
-- trusted reference), not values this implementation itself produced.

Testing.Assert_close(SpecialFunctions.Log_gamma(100), 359.1342053695754, 1e-9, "Log_gamma(100)")
Testing.Assert_close(SpecialFunctions.Log_gamma(170), 701.437263808737, 1e-9, "Log_gamma(170)")
Testing.Assert_close(SpecialFunctions.Log_gamma(0.3), 1.0957979948180752, 1e-9, "Log_gamma(0.3)")
Testing.Assert_close(SpecialFunctions.Log_gamma(-4.5), -2.8130840817693166, 1e-9, "Log_gamma(-4.5)")

-- Internal consistency: Log_gamma(x) == log(|Gamma(x)|) for moderate x
-- (where Gamma(x) itself doesn't overflow).
for _, x in ipairs({1.2, 2.5, 4.4, -0.5, -3.3, 7}) do
    Testing.Assert_close(SpecialFunctions.Log_gamma(x), math.log(math.abs(SpecialFunctions.Gamma(x))), 1e-8,
        "Log_gamma(x) == log(|Gamma(x)|) at x=" .. x)
end

Testing.Assert_error(function() return SpecialFunctions.Log_gamma(-2) end,
    "Log_gamma: x must not be a nonpositive integer")

-- ===== The overflow boundary: Gamma legitimately overflows to math.huge
-- once the true value exceeds a double's range (~Gamma(172) and up),
-- while Log_gamma -- staying in log-space the whole time -- does not. =====

assert(SpecialFunctions.Gamma(200) == math.huge, "Gamma(200) must overflow to +inf (the true value exceeds a double)")
Testing.Assert_close(SpecialFunctions.Log_gamma(200), 857.9336698258575, 1e-9,
    "Log_gamma(200) stays finite and accurate even where Gamma itself overflows")
