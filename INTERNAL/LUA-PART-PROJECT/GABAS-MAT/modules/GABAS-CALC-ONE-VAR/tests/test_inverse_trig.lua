local CalcOneVar = require("GABAS-CALC-ONE-VAR.GABAS_CALC_ONE_VAR")
local Complex = require("GABAS-LINAL.modules.complex")
local Testing = require("GABAS-LINAL.modules.testing")

local D = CalcOneVar.Derivative_at_point

-- Real-domain values, unchanged from math.asin/acos/atan.
local f_asin = CalcOneVar.Compile_expression("asin(x)")
local f_acos = CalcOneVar.Compile_expression("acos(x)")
local f_atan = CalcOneVar.Compile_expression("atan(x)")
Testing.Assert_close(f_asin(0.5), math.asin(0.5), 1e-12, "asin(0.5)")
Testing.Assert_close(f_acos(0.5), math.acos(0.5), 1e-12, "acos(0.5)")
Testing.Assert_close(f_atan(1), math.pi / 4, 1e-12, "atan(1)")

-- Real input outside [-1,1] auto-promotes to Complex (asin/acos), instead
-- of silently returning nan.
assert(Complex.Is_complex(f_asin(2)), "asin(2) should promote to Complex")
assert(Complex.Is_complex(f_acos(2)), "acos(2) should promote to Complex")

-- Genuinely Complex argument.
local f_asin_j = CalcOneVar.Compile_expression("asin(j*x)")
assert(Complex.Is_complex(f_asin_j(1)), "asin(j) should be Complex")

-- Derivative rules (Dual-aware): d/dx asin(x) = 1/sqrt(1-x^2),
-- d/dx acos(x) = -1/sqrt(1-x^2), d/dx atan(x) = 1/(1+x^2).
Testing.Assert_close((D(f_asin, 0.5)), 1 / math.sqrt(1 - 0.5 ^ 2), 1e-9, "D[asin(x)] at 0.5")
Testing.Assert_close((D(f_acos, 0.5)), -1 / math.sqrt(1 - 0.5 ^ 2), 1e-9, "D[acos(x)] at 0.5")
Testing.Assert_close((D(f_atan, 1)), 1 / (1 + 1 ^ 2), 1e-9, "D[atan(x)] at 1")

-- Chain rule through a composition: d/dx[asin(x^2)] = 2x/sqrt(1-x^4).
local f_comp = CalcOneVar.Compile_expression("asin(x^2)")
local x0 = 0.3
Testing.Assert_close((D(f_comp, x0)), 2 * x0 / math.sqrt(1 - x0 ^ 4), 1e-9, "D[asin(x^2)] at 0.3")
