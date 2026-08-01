local CalcOneVar = require("GABAS-CALC-ONE-VAR.GABAS_CALC_ONE_VAR")
local Testing = require("GABAS-LINAL.modules.testing")

local PV = CalcOneVar.Integral_principal_value

-- Symmetric case: PV(1/x) on [-1,1] = 0 -- the folded phi(t) = 1/t + 1/(-t)
-- cancels to exactly 0 (technical document, section 33).
local v1, e1, s1 = PV(function(x) return 1 / x end, -1, 1, 0)
Testing.Assert_close(v1, 0, 1e-9, "PV(1/x) on [-1,1]")
assert(s1 == "SUCCESS")

-- Asymmetric case: PV(1/x) on [-1,2] = ln(2) -- the symmetric part still
-- folds to 0, leaving only the ordinary leftover tail [1,2].
local v2, e2, s2 = PV(function(x) return 1 / x end, -1, 2, 0)
Testing.Assert_close(v2, math.log(2), 1e-6, "PV(1/x) on [-1,2]")
assert(s2 == "SUCCESS")

-- A regular part combined with a singular part: PV((x+1)/x) on [-1,1] = 2.
local v3, e3, s3 = PV(function(x) return (x + 1) / x end, -1, 1, 0)
Testing.Assert_close(v3, 2, 1e-6, "PV((x+1)/x) on [-1,1]")
assert(s3 == "SUCCESS")

-- Singularity not at the origin.
local v4, e4, s4 = PV(function(x) return 1 / (x - 1) end, 0, 2, 1)
Testing.Assert_close(v4, 0, 1e-6, "PV(1/(x-1)) on [0,2]")
assert(s4 == "SUCCESS")

-- An EVEN-type singularity (1/x^2) has no genuine principal value -- the
-- folding trick cannot cancel it (phi(t) = 2/t^2, still divergent), and
-- Integral_definite must not report a false SUCCESS for it.
local _, _, s5 = PV(function(x) return 1 / (x * x) end, -1, 1, 0)
assert(s5 ~= "SUCCESS", "1/x^2 has no principal value; must not report SUCCESS")

-- Input validation: c must lie strictly between a and b.
Testing.Assert_error(function() return PV(function(x) return 1 / x end, -1, 1, 1) end,
    "c must lie strictly between a and b")
Testing.Assert_error(function() return PV(function(x) return 1 / x end, -1, 1, -2) end,
    "c must lie strictly between a and b")
Testing.Assert_error(function() return PV(42, -1, 1, 0) end, "Integral_principal_value: f")
