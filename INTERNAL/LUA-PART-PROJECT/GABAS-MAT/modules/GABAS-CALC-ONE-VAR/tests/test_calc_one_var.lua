local CalcOneVar = require("GABAS-CALC-ONE-VAR.GABAS_CALC_ONE_VAR")
local Complex = require("GABAS-LINAL.modules.complex")
local Testing = require("GABAS-LINAL.modules.testing")

local D = CalcOneVar.Derivative_at_point
local f = CalcOneVar.Compile_expression("x^2")

-- Derivative_at_point's x0 must be exclusively a real number -- either of
-- Lua's own integer or float subtypes (math.type), signed or zero -- and
-- reject everything else: strings, booleans, nil, tables, NaN, +-infinity,
-- and Complex values. Enforced by Core.assert_finite_number (see
-- modules/core.lua), which is what type(x) == "number" already covers in
-- Lua 5.4 regardless of the integer/float subtype.

-- Accepted: integer or float, positive, negative, and zero.
Testing.Assert_close(D(f, 5), 10, 1e-9, "D(f,5)")
Testing.Assert_close(D(f, -5), -10, 1e-9, "D(f,-5)")
Testing.Assert_close(D(f, 3.14), 6.28, 1e-9, "D(f,3.14)")
Testing.Assert_close(D(f, -3.14), -6.28, 1e-9, "D(f,-3.14)")
Testing.Assert_close(D(f, 0), 0, 1e-9, "D(f,0)")
Testing.Assert_close(D(f, 0.0), 0, 1e-9, "D(f,0.0)")
Testing.Assert_close(D(f, -1e-10), -2e-10, 1e-9, "D(f,-1e-10)")

-- Rejected: anything that isn't a real, finite number.
Testing.Assert_error(function() return D(f, "5") end, "Derivative_at_point: x0")
Testing.Assert_error(function() return D(f, true) end, "Derivative_at_point: x0")
Testing.Assert_error(function() return D(f, nil) end, "Derivative_at_point: x0")
Testing.Assert_error(function() return D(f, {}) end, "Derivative_at_point: x0")
Testing.Assert_error(function() return D(f, 0 / 0) end, "Derivative_at_point: x0")
Testing.Assert_error(function() return D(f, math.huge) end, "Derivative_at_point: x0")
Testing.Assert_error(function() return D(f, -math.huge) end, "Derivative_at_point: x0")
Testing.Assert_error(function() return D(f, Complex.Complex(1, 2)) end, "Derivative_at_point: x0")
