local CalcOneVar = require("GABAS-CALC-ONE-VAR.GABAS_CALC_ONE_VAR")
local Testing = require("GABAS-LINAL.modules.testing")

-- f(x) = x^3 - 3x on [-3,3]: f'(x) = 3x^2-3 = 3(x-1)(x+1), critical
-- points at x=-1 (local MAXIMUM, f(-1)=2) and x=1 (local MINIMUM,
-- f(1)=-2) -- verified by hand, not just by this implementation.
-- Increasing on (-3,-1) and (1,3); decreasing on (-1,1).
local CUBIC = "x**3 - 3*x"

-- ===== Critical_points =====

local points, status = CalcOneVar.Critical_points(CUBIC, -3, 3)
assert(status == "SUCCESS")
assert(#points == 2, "x^3-3x on [-3,3] has exactly 2 critical points")
Testing.Assert_close(points[1].x, -1, 1e-6, "first critical point at x=-1")
assert(points[1].kind == "MAXIMUM", "x=-1 is a local maximum")
Testing.Assert_close(points[1].f_x, 2, 1e-6, "f(-1) = 2")
Testing.Assert_close(points[2].x, 1, 1e-6, "second critical point at x=1")
assert(points[2].kind == "MINIMUM", "x=1 is a local minimum")
Testing.Assert_close(points[2].f_x, -2, 1e-6, "f(1) = -2")

Testing.Assert_error(function() return CalcOneVar.Critical_points(42, -3, 3) end,
    "Critical_points: f")
Testing.Assert_error(function() return CalcOneVar.Critical_points(CUBIC, 3, -3) end,
    "Critical_points: a must be less than b")

-- ===== Local_maxima / Local_minima =====

local maxima = CalcOneVar.Local_maxima(CUBIC, -3, 3)
assert(#maxima == 1)
Testing.Assert_close(maxima[1].x, -1, 1e-6, "Local_maxima finds x=-1")

local minima = CalcOneVar.Local_minima(CUBIC, -3, 3)
assert(#minima == 1)
Testing.Assert_close(minima[1].x, 1, 1e-6, "Local_minima finds x=1")

-- ===== Increasing_intervals / Decreasing_intervals =====

local inc = CalcOneVar.Increasing_intervals(CUBIC, -3, 3)
assert(#inc == 2, "x^3-3x on [-3,3] increases on two separate pieces")
Testing.Assert_close(inc[1][1], -3, 1e-9, "first increasing interval starts at -3")
Testing.Assert_close(inc[1][2], -1, 1e-6, "first increasing interval ends at -1")
Testing.Assert_close(inc[2][1], 1, 1e-6, "second increasing interval starts at 1")
Testing.Assert_close(inc[2][2], 3, 1e-9, "second increasing interval ends at 3")

local dec = CalcOneVar.Decreasing_intervals(CUBIC, -3, 3)
assert(#dec == 1, "x^3-3x on [-3,3] decreases on exactly one piece")
Testing.Assert_close(dec[1][1], -1, 1e-6, "decreasing interval starts at -1")
Testing.Assert_close(dec[1][2], 1, 1e-6, "decreasing interval ends at 1")

Testing.Assert_error(function() return CalcOneVar.Increasing_intervals(42, -3, 3) end,
    "Increasing_intervals: f")

-- ===== Tangent_line =====

local slope, intercept, L = CalcOneVar.Tangent_line("x**2", 3)
Testing.Assert_close(slope, 6, 1e-9, "tangent slope of x^2 at x=3 is f'(3)=6")
Testing.Assert_close(intercept, -9, 1e-9, "tangent intercept")
Testing.Assert_close(L(3), 9, 1e-9, "tangent line matches f at the point of tangency")
Testing.Assert_close(L(0), -9, 1e-9, "tangent line evaluated elsewhere")

Testing.Assert_error(function() return CalcOneVar.Tangent_line(42, 3) end, "Tangent_line: f")

-- ===== Differential_approximation =====

local approx1, exact1, err1 = CalcOneVar.Differential_approximation("x**2", 3, 0.1)
Testing.Assert_close(approx1, 9.6, 1e-9, "linear approximation of (3.1)^2")
Testing.Assert_close(exact1, 9.61, 1e-9, "exact (3.1)^2")
Testing.Assert_close(err1, -0.01, 1e-9, "approximation error for a small delta_x")

-- A larger delta_x makes the linear approximation visibly worse --
-- confirms the error is a genuine measure of approximation quality, not
-- always trivially near 0.
local approx2, exact2, err2 = CalcOneVar.Differential_approximation("x**2", 3, 1)
Testing.Assert_close(approx2, 15, 1e-9, "linear approximation of 4^2 from x0=3")
Testing.Assert_close(exact2, 16, 1e-9, "exact 4^2")
Testing.Assert_close(err2, -1, 1e-9, "approximation error grows with delta_x")

Testing.Assert_error(function() return CalcOneVar.Differential_approximation(42, 3, 0.1) end,
    "Differential_approximation: f")
