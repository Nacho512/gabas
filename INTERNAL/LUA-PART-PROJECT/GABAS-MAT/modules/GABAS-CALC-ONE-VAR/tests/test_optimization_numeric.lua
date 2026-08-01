local CalcOneVar = require("GABAS-CALC-ONE-VAR.GABAS_CALC_ONE_VAR")
local Testing = require("GABAS-LINAL.modules.testing")

local GSMin = CalcOneVar.Golden_section_minimize
local GSMax = CalcOneVar.Golden_section_maximize
local LocalMin = CalcOneVar.Local_minimum
local LocalMax = CalcOneVar.Local_maximum
local GlobalMin = CalcOneVar.Global_minimum
local GlobalMax = CalcOneVar.Global_maximum

-- ===== Golden_section_minimize / Golden_section_maximize (derivative-
-- free, raw Lua functions with math.* are fine here) =====

local x1, f1, s1 = GSMin(function(x) return (x - 2) * (x - 2) end, 0, 5)
Testing.Assert_close(x1, 2, 1e-8, "Golden_section_minimize: (x-2)^2 on [0,5]")
Testing.Assert_close(f1, 0, 1e-8, "Golden_section_minimize: minimum value")
assert(s1 == "SUCCESS")

local x2, f2, s2 = GSMin(math.cos, 0, 2 * math.pi)
Testing.Assert_close(x2, math.pi, 1e-6, "Golden_section_minimize: cos(x) on [0,2pi]")
Testing.Assert_close(f2, -1, 1e-8, "Golden_section_minimize: cos minimum value")
assert(s2 == "SUCCESS")

-- Claude: deliberately NO large additive constant here (e.g. "+3") --
-- near the true optimum, a large constant swamps the tiny quadratic
-- term in floating-point subtraction (the quadratic term drops below
-- the constant's own ULP), making f numerically FLAT for a whole
-- neighborhood of x and starving golden-section's f(c) vs f(d)
-- comparisons of any real signal to work with. This is a genuine
-- floating-point precision limit of any comparison-based optimizer
-- measured this way, not specific to this implementation -- verified
-- directly: adding "+3" here alone widens the achievable x-precision
-- from ~1e-12 to ~1e-8.
local x3, f3, s3 = GSMax(function(x) return -(x - 2) * (x - 2) end, 0, 5)
Testing.Assert_close(x3, 2, 1e-8, "Golden_section_maximize: -(x-2)^2 on [0,5]")
Testing.Assert_close(f3, 0, 1e-8, "Golden_section_maximize: maximum value")
assert(s3 == "SUCCESS")

-- Accepts a symbolic-expression string.
local x4 = GSMin("(x - 2)**2", 0, 5)
Testing.Assert_close(x4, 2, 1e-8, "Golden_section_minimize with a string")

-- Input validation.
Testing.Assert_error(function() return GSMin(42, 0, 1) end, "Golden_section_minimize: f")
Testing.Assert_error(function() return GSMin(math.sin, 1, 1) end,
    "Golden_section_minimize: a and b must be distinct")

-- ===== Local_minimum / Local_maximum (derivative-based -- f must be
-- Dual-compatible, so a symbolic-expression string, exactly like Newton
-- in roots_numeric.lua) =====

local m1, mf1, ms1, md1 = LocalMin("x**2", 1)
Testing.Assert_close(m1, 0, 1e-8, "Local_minimum: x^2 from x0=1")
Testing.Assert_close(mf1, 0, 1e-8, "Local_minimum: minimum value")
assert(ms1 == "SUCCESS")
assert(md1 ~= nil)

local m2, mf2, ms2 = LocalMax("-(x**2)", 1)
Testing.Assert_close(m2, 0, 1e-8, "Local_maximum: -x^2 from x0=1")
Testing.Assert_close(mf2, 0, 1e-8, "Local_maximum: maximum value")
assert(ms2 == "SUCCESS")

-- x^3 at x0 near 0: f'(x)=3x^2=0 at x=0, but it is an INFLECTION point,
-- not a genuine extremum (f' does not change sign -- it is >= 0 on both
-- sides) -- must be reported honestly, never silently mislabeled.
local _, _, ms3 = LocalMin("x**3", 0.5)
assert(ms3 == "NOT_A_MINIMUM", "x=0 on x^3 is an inflection point, not a minimum")
local _, _, ms4 = LocalMax("x**3", 0.5)
assert(ms4 == "NOT_A_MAXIMUM", "x=0 on x^3 is an inflection point, not a maximum")

-- Input validation.
Testing.Assert_error(function() return LocalMin(42, 1) end, "Local_minimum: f")
Testing.Assert_error(function() return LocalMax(42, 1) end, "Local_maximum: f")

-- ===== Global_minimum / Global_maximum (derivative-based, same Dual-
-- compatibility requirement) =====
--
-- f(x) = x^4 - 4x^2 on [-3,3] -- NOT unimodal: a local max at x=0
-- (f(0)=0) and two symmetric local minima at x=+-sqrt(2) (f=-4), with
-- the interval endpoints even higher than the local max (f(+-3)=45).
-- This exercises both halves of the extreme-value-theorem logic: the
-- true global minimum is an INTERIOR critical point, while the true
-- global maximum is at an ENDPOINT, not the interior local max.

local gmin_x, gmin_f, gmin_s = GlobalMin("x**4 - 4*x**2", -3, 3)
Testing.Assert_close(gmin_f, -4, 1e-6, "Global_minimum: x^4-4x^2 on [-3,3]")
Testing.Assert_close(math.abs(gmin_x), math.sqrt(2), 1e-4, "Global_minimum: at +-sqrt(2)")
assert(gmin_s == "SUCCESS")

local gmax_x, gmax_f, gmax_s = GlobalMax("x**4 - 4*x**2", -3, 3)
Testing.Assert_close(gmax_f, 45, 1e-6, "Global_maximum: the true max is at an endpoint, not x=0")
Testing.Assert_close(math.abs(gmax_x), 3, 1e-9, "Global_maximum: at x=+-3")
assert(gmax_s == "SUCCESS")

-- A simple unimodal case, as a sanity cross-check against
-- Golden_section_minimize on the same problem.
local simple_x, simple_f = GlobalMin("(x - 2)**2", 0, 5)
Testing.Assert_close(simple_x, 2, 1e-4, "Global_minimum agrees with Golden_section_minimize")
Testing.Assert_close(simple_f, 0, 1e-6, "Global_minimum agrees with Golden_section_minimize (value)")

-- Input validation.
Testing.Assert_error(function() return GlobalMin(42, -3, 3) end, "Global_minimum: f")
Testing.Assert_error(function() return GlobalMin("x**2", 3, -3) end, "Global_minimum: a must be less than b")
