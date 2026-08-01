-- [Core]
-- Internal infrastructure for GABAS-DISCRETE-MATH's own submodules.
-- Nothing here is meant to be part of the public API -- every function
-- here is plumbing, not something an end user calls directly.
--
-- Claude: reuses GABAS-LINAL's own core.lua for every validator that
-- isn't specific to this domain (finite-number/integer checks, non-
-- blank-string checks, ...) rather than re-deriving them -- same pattern
-- GABAS-CALC-ONE-VAR's core.lua already established. What IS added here
-- is specific to this module: recognizing a valid numeral base (an
-- integer from 2 to 36, the largest base representable with the digits
-- 0-9 followed by the letters A-Z).
local LinalCore = require("GABAS-LINAL.modules.core")

local MIN_BASE, MAX_BASE = 2, 36

local function is_valid_base(base)
    return LinalCore.is_integer(base) and base >= MIN_BASE and base <= MAX_BASE
end

local function assert_valid_base(base, label)
    assert(is_valid_base(base),
        label .. " must be an integer base between " .. MIN_BASE .. " and " .. MAX_BASE .. ".")
end

return {
    -- Generic validators, reused from GABAS-LINAL's core rather than
    -- re-derived (same behavior, same messages policy).
    MACHINE_EPSILON = LinalCore.MACHINE_EPSILON,
    is_finite_number = LinalCore.is_finite_number,
    is_integer = LinalCore.is_integer,
    is_positive_number = LinalCore.is_positive_number,
    is_nonneg_number = LinalCore.is_nonneg_number,
    is_positive_integer = LinalCore.is_positive_integer,
    is_nonneg_integer = LinalCore.is_nonneg_integer,
    is_nonblank_string = LinalCore.is_nonblank_string,
    assert_finite_number = LinalCore.assert_finite_number,
    assert_integer = LinalCore.assert_integer,
    assert_positive_number = LinalCore.assert_positive_number,
    assert_nonneg_number = LinalCore.assert_nonneg_number,
    assert_positive_integer = LinalCore.assert_positive_integer,
    assert_nonneg_integer = LinalCore.assert_nonneg_integer,
    assert_nonblank_string = LinalCore.assert_nonblank_string,
    -- New here, specific to GABAS-DISCRETE-MATH.
    MIN_BASE = MIN_BASE,
    MAX_BASE = MAX_BASE,
    is_valid_base = is_valid_base,
    assert_valid_base = assert_valid_base,
}
