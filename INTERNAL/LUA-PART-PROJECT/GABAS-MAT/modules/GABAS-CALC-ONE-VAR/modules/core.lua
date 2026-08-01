-- [Core]
-- Internal infrastructure for GABAS-CALC-ONE-VAR's own submodules
-- (Compile_expression, Derivative_at_point, Root_any, Dual). Nothing here
-- is meant to be part of the public API -- every function here is
-- plumbing, not something an end user calls directly.
--
-- Claude: most of what user-input validation in THIS module needs is
-- already correct and already tested in GABAS-LINAL's own core.lua
-- (finite-number checks, integer checks, non-blank-string checks, ...) --
-- none of that is linear-algebra-specific, so this module requires it and
-- re-exports the generic pieces rather than re-deriving them. What IS
-- added here is specific to this domain: recognizing a Dual value (from
-- ./dual.lua, forward-mode automatic differentiation) as a legitimate
-- scalar, alongside plain reals and Complex values -- something
-- GABAS-LINAL's own core has no reason to know about.
--
-- Dependency direction matters here: this file may require Dual, but
-- Dual must NEVER require this file back (or GABAS_CALC_ONE_VAR.lua, or
-- anything built on top of it) -- that would be circular. Dual stays
-- self-contained, exactly like Complex.
local LinalCore = require("GABAS-LINAL.modules.core")
local Dual = require("GABAS-CALC-ONE-VAR.modules.dual")

-- True iff x is a finite real number, a finite Complex value, or a Dual
-- value whose own .value/.derivative are (recursively) finite scalars --
-- "is this one legitimate value that could flow through a compiled
-- expression", whether or not a derivative is being propagated.
local function is_finite_scalar_or_dual(x)
    if Dual.Is_dual(x) then
        return is_finite_scalar_or_dual(x.value) and is_finite_scalar_or_dual(x.derivative)
    end
    return LinalCore.is_finite_scalar(x)
end

local function assert_finite_scalar_or_dual(x, label)
    assert(is_finite_scalar_or_dual(x),
        label .. " must be a finite number, a finite Complex value, or a finite Dual value.")
end

-- Derivative_at_point/Verify_derivative_at_point both need "f is actually
-- callable" as their very first check, before ever trying to evaluate it.
local function is_callable(x)
    return type(x) == "function"
end

local function assert_callable(x, label)
    assert(is_callable(x), label .. " must be a function.")
end

return {
    -- Generic validators, reused from GABAS-LINAL's core rather than
    -- re-derived (same behavior, same messages policy).
    MACHINE_EPSILON = LinalCore.MACHINE_EPSILON,
    DEFAULT_ABS_TOL = LinalCore.DEFAULT_ABS_TOL,
    DEFAULT_REL_TOL = LinalCore.DEFAULT_REL_TOL,
    is_finite_number = LinalCore.is_finite_number,
    is_finite_scalar = LinalCore.is_finite_scalar,
    is_integer = LinalCore.is_integer,
    is_positive_number = LinalCore.is_positive_number,
    is_nonneg_number = LinalCore.is_nonneg_number,
    is_positive_integer = LinalCore.is_positive_integer,
    is_nonneg_integer = LinalCore.is_nonneg_integer,
    is_zero = LinalCore.is_zero,
    is_nonblank_string = LinalCore.is_nonblank_string,
    assert_finite_number = LinalCore.assert_finite_number,
    assert_finite_scalar = LinalCore.assert_finite_scalar,
    assert_integer = LinalCore.assert_integer,
    assert_positive_number = LinalCore.assert_positive_number,
    assert_nonneg_number = LinalCore.assert_nonneg_number,
    assert_positive_integer = LinalCore.assert_positive_integer,
    assert_nonneg_integer = LinalCore.assert_nonneg_integer,
    assert_nonblank_string = LinalCore.assert_nonblank_string,
    -- New here, specific to GABAS-CALC-ONE-VAR.
    is_finite_scalar_or_dual = is_finite_scalar_or_dual,
    assert_finite_scalar_or_dual = assert_finite_scalar_or_dual,
    is_callable = is_callable,
    assert_callable = assert_callable,
}
