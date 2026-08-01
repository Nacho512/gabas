local CalcOneVar = require("GABAS-CALC-ONE-VAR.GABAS_CALC_ONE_VAR")
local Testing = require("GABAS-LINAL.modules.testing")

local Osc = CalcOneVar.Integral_oscillatory

-- Filon-type quadrature (technical document, section 14): integral of
-- g(x)*cos(omega*x) or g(x)*sin(omega*x) over a finite interval, with
-- omega possibly large. Reference values below are independent
-- high-precision numerical integrals (mpmath, 40 digits), not values this
-- implementation itself produced.

-- x*cos(100x) on [0,1] -- a genuinely rapid oscillation (~16 periods)
-- against a simple linear amplitude.
local v1, e1, s1, d1 = Osc(function(x) return x end, 0, 1, 100, "cos")
Testing.Assert_close(v1, -0.005077424523868819543155382253202770067099, 1e-6,
    "x*cos(100x) on [0,1]")
assert(s1 == "SUCCESS")
assert(d1.levels ~= nil and d1.panels ~= nil)

-- sin(x)*sin(50x) on [0,pi] -- product-to-sum identity makes this
-- exactly 0 (cos(49x) and cos(51x) both integrate to 0 over a symmetric
-- number of half-periods).
local v2, e2, s2 = Osc(function(x) return math.sin(x) end, 0, math.pi, 50, "sin")
Testing.Assert_close(v2, 0, 1e-6, "sin(x)*sin(50x) on [0,pi]")
assert(s2 == "SUCCESS")

-- e^-x*sin(30x) on [0,pi] -- closed form via the standard
-- exponential-times-sine antiderivative.
local v3, e3, s3 = Osc(function(x) return math.exp(-x) end, 0, math.pi, 30, "sin")
Testing.Assert_close(v3, 0.03185747219987439789874302761914335145586, 1e-6,
    "e^-x*sin(30x) on [0,pi]")
assert(s3 == "SUCCESS")

-- cos(200x) on [0,2*pi] with a CONSTANT amplitude -- exactly 200 full
-- periods, so the integral is 0 to essentially machine precision.
local v4, e4, s4 = Osc(function(x) return 1 end, 0, 2 * math.pi, 200, "cos")
Testing.Assert_close(v4, 0, 1e-6, "cos(200x) on [0,2*pi], 200 full periods")
assert(s4 == "SUCCESS")

-- cos(17x) on [0.3,2.7] with a constant amplitude -- Filon is exact for a
-- constant amplitude (the quadratic interpolant of a constant IS the
-- constant), so this should converge extremely tightly.
local v5, e5, s5 = Osc(function(x) return 1 end, 0.3, 2.7, 17, "cos")
Testing.Assert_close(v5, 0.10977903795733741, 1e-9, "cos(17x) on [0.3,2.7], constant amplitude")
assert(s5 == "SUCCESS")

-- Cross-check against the ordinary adaptive engine (Integral_definite),
-- which can also resolve a moderately oscillatory integral if allowed
-- enough subintervals -- independent verification, not just agreement
-- with itself.
local v6 = Osc(function(x) return math.exp(x) end, 0, 1, 40, "cos")
local v6_ref = CalcOneVar.Integral_definite(function(x) return math.exp(x) * math.cos(40 * x) end,
    0, 1, {max_intervals = 2000})
Testing.Assert_close(v6, v6_ref, 1e-6, "exp(x)*cos(40x) on [0,1], cross-checked against Integral_definite")

-- Accepts a symbolic-expression string for g directly (Resolve_function),
-- same as every other Integral_* function in this module.
local v7 = Osc("x**2", 0, 1, 60, "sin")
local v7_direct = Osc(function(x) return x * x end, 0, 1, 60, "sin")
Testing.Assert_close(v7, v7_direct, 1e-12, "Integral_oscillatory with a string amplitude")

-- Input validation.
Testing.Assert_error(function() return Osc(42, 0, 1, 10, "cos") end, "Integral_oscillatory: g")
Testing.Assert_error(function() return Osc(math.sin, 1, 0, 10, "cos") end,
    "Integral_oscillatory: a must be less than b")
Testing.Assert_error(function() return Osc(math.sin, 0, 1, 0, "cos") end,
    "Integral_oscillatory: omega must be nonzero")
Testing.Assert_error(function() return Osc(math.sin, 0, 1, 10, "tan") end,
    "Integral_oscillatory: kind must be")
Testing.Assert_error(function() return Osc(math.sin, 0, 1, 10, "cos", {max_levels = -1}) end,
    "Integral_oscillatory: opts.max_levels")
