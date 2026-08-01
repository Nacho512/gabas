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

-- ===== Beta / Log_beta =====
--
-- Reference values below are mpmath's own beta() at 30-digit working
-- precision (an independent, trusted reference distinct from this
-- module's own Log_gamma-based formula), plus a cross-check against
-- Beta's DEFINING INTEGRAL for a few values -- two independent
-- verifications, not just one formula checked against itself.

Testing.Assert_close(SpecialFunctions.Beta(1, 1), 1, 1e-12, "Beta(1,1) = 1")
Testing.Assert_close(SpecialFunctions.Beta(0.5, 0.5), math.pi, 1e-9, "Beta(0.5,0.5) = pi")
Testing.Assert_close(SpecialFunctions.Beta(2, 3), 1 / 12, 1e-12, "Beta(2,3) = 1/12")
Testing.Assert_close(SpecialFunctions.Beta(10, 20), 4.9925087406346777e-09, 1e-9, "Beta(10,20)")
Testing.Assert_close(SpecialFunctions.Beta(2.5, 3.7), 0.032727368606257838, 1e-9, "Beta(2.5,3.7)")
Testing.Assert_close(SpecialFunctions.Beta(0.1, 0.1), 19.714639489050161, 1e-8, "Beta(0.1,0.1)")

-- The overflow case: Gamma(100)*Gamma(100)/Gamma(200) computed naively
-- would divide by +inf (Gamma(200) already overflows a double) and give
-- a silently WRONG 0 -- Beta(100,100) is actually an ordinary, tiny but
-- nonzero positive number. Verified this really is what the naive
-- formula does, in Python, before writing Log_beta this way.
Testing.Assert_close(SpecialFunctions.Beta(100, 100), 2.2087606931995026e-61, 1e-70, "Beta(100,100), the overflow case")

-- Cross-check against the defining integral, via a fine Riemann sum
-- (independent of Log_gamma entirely) -- not exact quadrature, so a
-- looser tolerance.
local function beta_by_integral(a, b, n)
    local sum, h = 0, 1 / n
    for i = 1, n - 1 do
        local t = i * h
        sum = sum + (t ^ (a - 1)) * ((1 - t) ^ (b - 1))
    end
    return sum * h
end
Testing.Assert_close(SpecialFunctions.Beta(2, 3), beta_by_integral(2, 3, 200000), 1e-4,
    "Beta(2,3) matches its own defining integral, computed independently")
Testing.Assert_close(SpecialFunctions.Beta(5, 5), beta_by_integral(5, 5, 200000), 1e-6,
    "Beta(5,5) matches its own defining integral, computed independently")

-- Log_beta internal consistency.
for _, ab in ipairs({{1, 1}, {2, 3}, {10, 20}, {100, 100}}) do
    Testing.Assert_close(SpecialFunctions.Log_beta(ab[1], ab[2]), math.log(SpecialFunctions.Beta(ab[1], ab[2])), 1e-8,
        "Log_beta == log(Beta) at a=" .. ab[1] .. ",b=" .. ab[2])
end

-- Input validation: a, b must both be strictly positive.
Testing.Assert_error(function() return SpecialFunctions.Beta(0, 1) end, "Beta: a must be positive")
Testing.Assert_error(function() return SpecialFunctions.Beta(1, -1) end, "Beta: b must be positive")
Testing.Assert_error(function() return SpecialFunctions.Beta("1", 1) end, "Beta: a")
Testing.Assert_error(function() return SpecialFunctions.Log_beta(0, 1) end, "Log_beta: a must be positive")
