local CalcOneVar = require("GABAS-CALC-ONE-VAR.GABAS_CALC_ONE_VAR")
local Complex = require("GABAS-LINAL.modules.complex")
local Testing = require("GABAS-LINAL.modules.testing")

local I = CalcOneVar.Integral_definite

-- Basic analytical integrals (technical document, section 27).
Testing.Assert_close((I(function(x) return x ^ 19 end, 0, 1)), 1 / 20, 1e-8, "x^19 on [0,1]")
Testing.Assert_close((I(math.exp, 0, 1)), math.exp(1) - 1, 1e-8, "exp(x) on [0,1]")
Testing.Assert_close((I(math.sin, 0, math.pi)), 2, 1e-8, "sin(x) on [0,pi]")
Testing.Assert_close((I(math.cos, 0, 2 * math.pi)), 0, 1e-7, "cos(x) on [0,2pi]")
Testing.Assert_close((I(function(x) return 1 / x end, 1, math.exp(1))), 1, 1e-8, "1/x on [1,e]")
Testing.Assert_close((I(function(x) return 1 / (1 + x * x) end, 0, 1)), math.pi / 4, 1e-8, "1/(1+x^2) on [0,1]")

-- Complex-valued integrand (section 28): integral of e^(ix) on [0,pi] = 2i.
local f_complex = CalcOneVar.Compile_expression("exp(j*x)")
local v_complex, err_complex, status_complex = I(f_complex, 0, math.pi)
assert(Complex.Is_complex(v_complex), "exp(j*x) integral should be Complex-valued")
Testing.Assert_close(v_complex.re, 0, 1e-8, "Re(integral of e^(ix))")
Testing.Assert_close(v_complex.im, 2, 1e-8, "Im(integral of e^(ix))")
assert(status_complex == "SUCCESS", "exp(j*x) integral should converge")

-- Bound handling (section 29): equal limits return 0 without evaluating
-- f (must not error even though f is undefined at that exact point).
local value_equal, error_equal, status_equal, diag_equal =
    I(function(x) return 1 / (x - 3) end, 3, 3)
assert(value_equal == 0 and error_equal == 0, "equal limits must give exactly 0")
assert(status_equal == "SUCCESS", "equal limits must report SUCCESS")
assert(diag_equal.evaluations == 0, "equal limits must not evaluate f at all")

-- Reversed limits: I(b,a) = -I(a,b).
Testing.Assert_close((I(math.sin, 0, math.pi)), -(I(math.sin, math.pi, 0)), 1e-8, "orientation reversal")

-- Metamorphic tests (section 30): linearity and interval additivity.
local v_lin = I(function(x) return 2 * math.sin(x) + 3 * math.cos(x) end, 0, math.pi)
local v_lin_expected = 2 * (I(math.sin, 0, math.pi)) + 3 * (I(math.cos, 0, math.pi))
Testing.Assert_close(v_lin, v_lin_expected, 1e-8, "linearity")

local v_full = I(math.exp, 0, 2)
local v_split = (I(math.exp, 0, 1)) + (I(math.exp, 1, 2))
Testing.Assert_close(v_full, v_split, 1e-8, "interval additivity")

-- Honest non-convergence: sin(x^7+e^x) on [2,10] is genuinely highly
-- oscillatory near the top of the range (the inner argument's derivative
-- grows from ~455 to ~7 million radians per unit x) -- exactly the
-- specialized-method case explicitly deferred in this first
-- implementation. The algorithm must report a non-SUCCESS status rather
-- than a false convergence; on a mild sub-range it must still succeed
-- cleanly.
local f_oscillatory = CalcOneVar.Compile_expression("sin(x**7+e(x))")
local _, _, status_mild = I(f_oscillatory, 2, 2.5)
assert(status_mild == "SUCCESS", "a mild sub-range should converge cleanly")
local _, _, status_hard = I(f_oscillatory, 2, 10, {max_intervals = 200})
assert(status_hard ~= "SUCCESS", "a highly oscillatory range must not report false convergence")

-- Input validation, mirroring Derivative_at_point's contract.
Testing.Assert_error(function() return I("not a function", 0, 1) end, "Integral_definite: f")
Testing.Assert_error(function() return I(math.sin, "0", 1) end, "Integral_definite: a")
Testing.Assert_error(function() return I(math.sin, 0, 1, {abs_tol = -1}) end, "Integral_definite: opts.abs_tol")
Testing.Assert_error(function() return I(math.sin, 0, 1, {abs_tol = 0, rel_tol = 0}) end,
    "at least one of opts.abs_tol/opts.rel_tol")
