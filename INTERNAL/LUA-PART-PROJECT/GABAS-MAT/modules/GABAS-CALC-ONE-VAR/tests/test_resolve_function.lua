local CalcOneVar = require("GABAS-CALC-ONE-VAR.GABAS_CALC_ONE_VAR")
local Testing = require("GABAS-LINAL.modules.testing")

-- Every public function that takes f accepts EITHER a symbolic-expression
-- string (compiled automatically via Compile_expression's Resolve_function)
-- OR an already-compiled function, so a caller never has to remember to
-- call Compile_expression first.

Testing.Assert_close((CalcOneVar.Derivative_at_point("sin(x**2)", 0.5)),
    math.cos(0.5 ^ 2) * 2 * 0.5, 1e-9, "Derivative_at_point with a string")

Testing.Assert_close((CalcOneVar.Integral_definite("sin(x)", 0, math.pi)), 2, 1e-8,
    "Integral_definite with a string")

Testing.Assert_close((CalcOneVar.Integral_infinite("exp(-x)", 0, math.huge)), 1, 1e-6,
    "Integral_infinite with a string")

Testing.Assert_close((CalcOneVar.Integral_principal_value("1/x", -1, 1, 0)), 0, 1e-9,
    "Integral_principal_value with a string")

Testing.Assert_close((CalcOneVar.Integral_tanh_sinh("x**(-0.5)", 0, 1)), 2, 1e-8,
    "Integral_tanh_sinh with a string")

-- Passing an already-compiled function still works unchanged (backward
-- compatible with every earlier example in this test suite).
local f = CalcOneVar.Compile_expression("x**2")
Testing.Assert_close((CalcOneVar.Derivative_at_point(f, 3)), 6, 1e-9,
    "Derivative_at_point with an already-compiled function")

-- A malformed expression string still fails to compile, with
-- Compile_expression's own error, not a generic "must be a function".
Testing.Assert_error(function() return CalcOneVar.Derivative_at_point("sin((x", 0.5) end,
    "Compile_Expression")

-- Anything that is neither a string nor a function is rejected.
Testing.Assert_error(function() return CalcOneVar.Derivative_at_point(42, 0.5) end,
    "Derivative_at_point: f")
Testing.Assert_error(function() return CalcOneVar.Derivative_at_point(true, 0.5) end,
    "Derivative_at_point: f")
