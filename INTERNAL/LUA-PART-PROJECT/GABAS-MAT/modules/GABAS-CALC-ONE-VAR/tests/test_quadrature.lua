local Quadrature = require("GABAS-CALC-ONE-VAR.modules.quadrature")
local Complex = require("GABAS-LINAL.modules.complex")
local Testing = require("GABAS-LINAL.modules.testing")

-- The 21-point Kronrod rule alone, exercised directly (no subdivision) --
-- verifies the node/weight tables themselves are correct, independent of
-- the global adaptive loop in GABAS_CALC_ONE_VAR.lua.
Testing.Assert_close((Quadrature.Local_Gauss_Kronrod_21(function(x) return 1 end, -1, 1)), 2, 1e-13,
    "integral of 1 on [-1,1]")
Testing.Assert_close((Quadrature.Local_Gauss_Kronrod_21(function(x) return x ^ 19 end, 0, 1)), 1 / 20, 1e-12,
    "x^19 on [0,1] (exact for the embedded 10-point Gauss rule)")
Testing.Assert_close((Quadrature.Local_Gauss_Kronrod_21(math.exp, 0, 1)), math.exp(1) - 1, 1e-13,
    "exp(x) on [0,1]")

-- Complex-valued integrand.
local j = Complex.Complex(0, 1)
local v = Quadrature.Local_Gauss_Kronrod_21(function(x) return Complex.Exp(j * x) end, 0, math.pi)
assert(Complex.Is_complex(v), "e^(ix) result should be Complex")
Testing.Assert_close(v.re, 0, 1e-8, "Re(e^(ix)) integral")
Testing.Assert_close(v.im, 2, 1e-8, "Im(e^(ix)) integral")

-- Node-evaluation failure is reported, not silently propagated as NaN.
local result, err, resabs, resasc, evals, message =
    Quadrature.Local_Gauss_Kronrod_21(function(x) return 1 / x end, -1, 1)
assert(result == nil, "1/x on [-1,1] (pole exactly at the center node) must fail, not return a number")
assert(type(message) == "string" and message:find("infinite", 1, true), "failure message must mention the cause")
