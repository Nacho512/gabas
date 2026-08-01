local CalcOneVar = require("GABAS-CALC-ONE-VAR.GABAS_CALC_ONE_VAR")
local Testing = require("GABAS-LINAL.modules.testing")

local Bisection = CalcOneVar.Bisection
local Newton = CalcOneVar.Newton
local Secant = CalcOneVar.Secant

-- Reference roots below are independent high-precision values (mpmath,
-- 30 digits), not values these implementations themselves produced.
local SQRT2 = 1.41421356237309504880168872421
local CUBIC_ROOT = 1.52137970680456756960408083225 -- root of x^3 - x - 2
local DOTTIE = 0.739085133215160641655312087674 -- root of cos(x) - x
local LN3 = 1.09861228866810969139524523692

-- ===== Bisection =====

local r1, e1, s1, d1 = Bisection(function(x) return x * x - 2 end, 1, 2)
Testing.Assert_close(r1, SQRT2, 1e-10, "Bisection: x^2-2 on [1,2]")
assert(s1 == "SUCCESS")
assert(e1 <= 1e-10)

local r2, e2, s2 = Bisection(function(x) return math.cos(x) - x end, 0, 1)
Testing.Assert_close(r2, DOTTIE, 1e-10, "Bisection: cos(x)-x on [0,1]")
assert(s2 == "SUCCESS")

-- Root exactly at an endpoint -- detected immediately, 0 iterations.
local r3, e3, s3, d3 = Bisection(function(x) return x - 1 end, -1, 1)
Testing.Assert_close(r3, 1, 1e-15, "Bisection: exact root at endpoint b")
assert(s3 == "SUCCESS")
assert(d3.iterations == 0)

-- Reversed bracket (b < a) is still accepted.
local r4, e4, s4 = Bisection(function(x) return x * x - 2 end, 2, 1)
Testing.Assert_close(r4, SQRT2, 1e-10, "Bisection: reversed bracket")
assert(s4 == "SUCCESS")

-- Accepts a symbolic-expression string.
local r5 = Bisection("x**2 - 2", 1, 2)
Testing.Assert_close(r5, SQRT2, 1e-10, "Bisection with a string")

-- Input validation.
Testing.Assert_error(function() return Bisection(42, 1, 2) end, "Bisection: f")
Testing.Assert_error(function() return Bisection(function(x) return x end, 1, 1) end,
    "Bisection: a and b must be distinct")
Testing.Assert_error(function() return Bisection(function(x) return x * x + 1 end, -1, 1) end,
    "Bisection: f(a) and f(b) must have opposite signs")

-- ===== Newton =====
--
-- Newton seeds its argument with a Dual value (forward-mode automatic
-- differentiation), so every f passed to it must be built from the
-- library's own Dual-aware primitives -- a symbolic-expression string
-- (compiled via Compile_expression, which every Calc_one_var_* primitive
-- already dispatches on Dual for) rather than a raw Lua function calling
-- math.* directly (math.exp etc. only accept plain numbers, not Dual
-- values).

local n1, ne1, ns1, nd1 = Newton("x**2 - 2", 1)
Testing.Assert_close(n1, SQRT2, 1e-12, "Newton: x^2-2 from x0=1")
assert(ns1 == "SUCCESS")
assert(nd1.iterations <= 10, "Newton should converge quadratically, in a handful of iterations")

local n2, ne2, ns2 = Newton("x**3 - x - 2", 1.5)
Testing.Assert_close(n2, CUBIC_ROOT, 1e-10, "Newton: x^3-x-2 from x0=1.5")
assert(ns2 == "SUCCESS")

local n3, ne3, ns3 = Newton("exp(x) - 3", 1)
Testing.Assert_close(n3, LN3, 1e-10, "Newton: exp(x)-3 from x0=1")
assert(ns3 == "SUCCESS")

-- x0 is already an exact root, and f'(x0) is ALSO 0 there (f(x)=x^2 at
-- x0=0) -- must succeed immediately, not trip the min_derivative guard.
local n4, ne4, ns4, nd4 = Newton("x**2", 0)
Testing.Assert_close(n4, 0, 1e-15, "Newton: x0 already an exact root with f'(x0)=0 too")
assert(ns4 == "SUCCESS")
assert(nd4.iterations == 0)

-- f'(x0) is 0 and f(x0) is NOT 0 (x^2+5 has no real root) -- must report
-- DERIVATIVE_TOO_SMALL honestly, not silently divide by (near) zero.
local _, _, ns5, nd5 = Newton("x**2 + 5", 0)
assert(ns5 == "DERIVATIVE_TOO_SMALL")
Testing.Assert_close(nd5.residual, 5, 1e-15, "Newton: residual reported for the no-real-root case")

-- Classical period-2 cycle (x^3-2x+2 from x0=0: 0 -> 1 -> 0 -> 1 -> ...,
-- verified by hand) -- Newton must never falsely report SUCCESS here.
local _, _, ns6 = Newton("x**3 - 2*x + 2", 0, {max_iterations = 20})
assert(ns6 ~= "SUCCESS", "Newton must not falsely converge on a genuine cycling case")

-- Accepts a symbolic-expression string.
local n7 = Newton("x**2 - 2", 1)
Testing.Assert_close(n7, SQRT2, 1e-10, "Newton with a string")

-- Input validation.
Testing.Assert_error(function() return Newton(42, 1) end, "Newton: f")
Testing.Assert_error(function() return Newton(function(x) return x end, 1, {tol = -1}) end,
    "Newton: opts.tol")

-- ===== Secant =====

local c1, ce1, cs1 = Secant(function(x) return x * x - 2 end, 1, 2)
Testing.Assert_close(c1, SQRT2, 1e-10, "Secant: x^2-2 from x0=1,x1=2")
assert(cs1 == "SUCCESS")

local c2, ce2, cs2 = Secant(function(x) return math.cos(x) - x end, 0, 1)
Testing.Assert_close(c2, DOTTIE, 1e-10, "Secant: cos(x)-x from x0=0,x1=1")
assert(cs2 == "SUCCESS")

-- Two starting points with EXACTLY equal f-values ((x-1)^2 is symmetric
-- about x=1, and 0/2 are equidistant from it) -- must report STALLED
-- rather than divide by zero.
local _, _, cs3, cd3 = Secant(function(x) return (x - 1) * (x - 1) end, 0, 2)
assert(cs3 == "STALLED")
assert(cd3.iterations == 0)

-- Accepts a symbolic-expression string.
local c4 = Secant("x**2 - 2", 1, 2)
Testing.Assert_close(c4, SQRT2, 1e-10, "Secant with a string")

-- Input validation.
Testing.Assert_error(function() return Secant(42, 1, 2) end, "Secant: f")
Testing.Assert_error(function() return Secant(function(x) return x end, 1, 1) end,
    "Secant: x0 and x1 must be distinct")
