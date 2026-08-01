local CalcOneVar = require("GABAS-CALC-ONE-VAR.GABAS_CALC_ONE_VAR")
local Testing = require("GABAS-LINAL.modules.testing")

local TS = CalcOneVar.Integral_tanh_sinh

-- Endpoint singularities (technical document, section 31) -- these
-- converge to machine precision.
local v1, e1, s1 = TS(function(x) return x ^ (-0.5) end, 0, 1)
Testing.Assert_close(v1, 2, 1e-9, "x^-0.5 on [0,1]")
assert(s1 == "SUCCESS")

local v2, e2, s2 = TS(math.log, 0, 1)
Testing.Assert_close(v2, -1, 1e-9, "ln(x) on [0,1]")
assert(s2 == "SUCCESS")

-- x^alpha near the alpha=-1 (non-integrable) boundary: still converges
-- comfortably at alpha=-0.9.
local v3, e3, s3 = TS(function(x) return x ^ (-0.9) end, 0, 1)
Testing.Assert_close(v3, 10, 1e-6, "x^-0.9 on [0,1]")
assert(s3 == "SUCCESS")

-- An INTERIOR singularity (here, |x|^-0.5 at x=0, the middle of
-- [-1,1]) is a different problem this function does not detect on its
-- own -- split at the singular point first, then sum both halves.
local left = TS(function(x) return math.abs(x) ^ (-0.5) end, -1, 0)
local right = TS(function(x) return math.abs(x) ^ (-0.5) end, 0, 1)
Testing.Assert_close(left + right, 4, 1e-6, "|x|^-0.5 on [-1,1], split at the singularity")

-- Genuinely divergent (alpha=-1, the harmonic-like case): must not
-- report SUCCESS, no matter how many levels are allowed.
local _, _, s4 = TS(function(x) return 1 / x end, 0, 1, {max_levels = 16})
assert(s4 ~= "SUCCESS", "the divergent integral of 1/x on [0,1] must not report SUCCESS")

-- A singularity close enough to the alpha=-1 boundary hits a genuine
-- precision limit that more levels do not fix -- must still report
-- honestly rather than a false SUCCESS.
local _, _, s5 = TS(function(x) return x ^ (-0.99) end, 0, 1, {max_levels = 20})
assert(s5 ~= "SUCCESS", "alpha=-0.99 should not falsely report SUCCESS at this precision")

-- An ordinary, non-singular integrand still converges correctly (and
-- quickly).
local v6, e6, s6, d6 = TS(math.sin, 0, math.pi)
Testing.Assert_close(v6, 2, 1e-9, "sin(x) on [0,pi] (no singularity)")
assert(s6 == "SUCCESS")
assert(d6.levels <= 8, "an ordinary integrand should converge in a handful of levels")

-- Input validation.
Testing.Assert_error(function() return TS(42, 0, 1) end, "Integral_tanh_sinh: f")
Testing.Assert_error(function() return TS(math.sin, 1, 0) end, "Integral_tanh_sinh: a must be less than b")

-- The interval's own center is checked strictly: a genuinely broken f
-- (not merely singular at an endpoint) must fail loudly, not silently
-- return a wrong number.
Testing.Assert_error(function()
    return TS(function() error("deliberately broken") end, 0, 1)
end, "the interval's center")
