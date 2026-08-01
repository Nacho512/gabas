-- [Dif Numeric]
-- Derivative_at_point uses forward-mode automatic differentiation via
-- Dual numbers (modules/dual.lua), not finite differences -- exact up
-- to floating-point rounding, no step-size parameter, no truncation
-- error. Every Calc_one_var_* primitive in compile_expression.lua
-- already dispatches on whether its argument is Dual, the same way it
-- already dispatches on Complex for "j" -- so differentiating a compiled
-- expression is just calling that SAME function with a Dual seed instead
-- of a plain number, and the derivative propagates through the existing
-- bytecode via operator overloading. Deliberately scoped to a REAL x only
-- (see dual.lua's header) -- differentiating with respect to a complex
-- variable itself (Jacobian/Wirtinger modes) is a different problem,
-- planned as its own separate module.
local Complex = require("GABAS-LINAL.modules.complex")
local Dual = require("GABAS-CALC-ONE-VAR.modules.dual")
local Core = require("GABAS-CALC-ONE-VAR.modules.core")

-- Claude: NaN-safe on both plain reals and Complex values -- used to
-- reject a derivative that silently went bad (e.g. from a domain issue an
-- individual primitive rule didn't explicitly guard) rather than handing
-- the caller a NaN with no explanation.
local function has_nan(v)
    if Complex.Is_complex(v) then
        return v.re ~= v.re or v.im ~= v.im
    end
    return type(v) == "number" and v ~= v
end

-- Derivative_at_point(f, x0): the numerical derivative of a real-to-real-
-- or-complex function f at a real point x0, via forward-mode automatic
-- differentiation (dual numbers) -- exact up to floating-point rounding,
-- with no step-size parameter and no finite-difference truncation error,
-- unlike ordinary numerical differentiation.
--
-- f is normally the function returned by Compile_expression, though any
-- Lua function of one argument works: the derivative comes from calling
-- that SAME f with a Dual seed instead of a plain number. Every
-- arithmetic operation and every Calc_one_var_* primitive f's body might
-- use already knows how to propagate a Dual value -- exactly the same way
-- they already know how to propagate a Complex one for "j" -- so no
-- separate parser or expression tree is needed; the derivative rides
-- through the SAME compiled bytecode Compile_expression already produced.
--
-- Deliberately scoped to a REAL x0 only (see the Dual module header) --
-- differentiating with respect to a complex variable itself is a
-- different problem (holomorphicity, branch cuts, Jacobian/Wirtinger
-- modes), planned as a separate module, not folded in here.
--
-- Returns derivative, value: derivative first, since that is this
-- function's whole purpose; value = f(x0) comes back too, already
-- computed along the way at no extra cost, useful e.g. to sanity-check
-- against a known value.
local function Derivative_at_point(f, x0)
    Core.assert_callable(f, "Derivative_at_point: f")
    Core.assert_finite_number(x0, "Derivative_at_point: x0")
    local seed = Dual.Dual(x0, 1)
    local ok, result = pcall(f, seed)
    assert(ok, "Derivative_at_point: error while evaluating f at x0 -- " .. tostring(result))
    if Dual.Is_dual(result) then
        local value, derivative = result.value, result.derivative
        assert(not has_nan(value), "Derivative_at_point: the expression evaluated to NaN at x0.")
        assert(not has_nan(derivative), "Derivative_at_point: the derivative evaluated to NaN at x0.")
        return derivative, value
    end
    -- f never actually used x (e.g. a constant expression like "7") --
    -- its derivative is trivially 0 everywhere.
    return 0, result
end

local function magnitude(v)
    if Complex.Is_complex(v) then return v:abs() end
    return math.abs(v)
end

-- Claude: an OPTIONAL, independent cross-check against a central finite
-- difference (technical document, section 15.1) -- NOT the primary
-- method (Derivative_at_point's dual-number result is already exact up
-- to rounding), just a cheap second opinion that would catch a mistake
-- in one of the hand-coded derivative rules above. h follows the standard
-- balance between central-difference truncation error (~h^2) and
-- floating-point round-off (~epsilon/h): h = cbrt(machine epsilon),
-- scaled by max(1, |x0|) so it stays meaningful for a large x0 too.
--
-- Returns agrees (boolean), finite_difference (the independent estimate).
local function Verify_derivative_at_point(f, x0, derivative, tolerance)
    Core.assert_callable(f, "Verify_derivative_at_point: f")
    Core.assert_finite_number(x0, "Verify_derivative_at_point: x0")
    Core.assert_finite_scalar_or_dual(derivative, "Verify_derivative_at_point: derivative")
    tolerance = tolerance or 1e-4
    Core.assert_nonneg_number(tolerance, "Verify_derivative_at_point: tolerance")
    local h = (Core.MACHINE_EPSILON ^ (1 / 3)) * math.max(1, math.abs(x0))
    local finite_difference = (f(x0 + h) - f(x0 - h)) / (2 * h)
    local agrees = magnitude(finite_difference - derivative) <= tolerance
    return agrees, finite_difference
end

return {
    Derivative_at_point = Derivative_at_point,
    Verify_derivative_at_point = Verify_derivative_at_point,
}
