-- [Core]
-- Internal infrastructure for GABAS-COMBINATORICS's own submodules.
-- Nothing here is meant to be part of the public API -- every function
-- here is plumbing, not something an end user calls directly.
--
-- Claude: reuses GABAS-LINAL's own core.lua for every validator this
-- module needs (finite-number/integer checks, ...) rather than re-
-- deriving them -- same pattern GABAS-CALC-ONE-VAR's and
-- GABAS-DISCRETE-MATH's own core.lua already established. Required
-- directly from GABAS-LINAL (not through GABAS-DISCRETE-MATH), since
-- combinatorics has no actual need for base-conversion-specific
-- validators -- each domain's core.lua depends on the true shared
-- foundation directly, keeping GABAS-DISCRETE-MATH and
-- GABAS-COMBINATORICS independent siblings rather than needlessly
-- coupled to each other.
--
-- Nothing domain-specific has been needed here yet -- Factorial,
-- Binomial_coefficient, Pascal_row/Pascal_triangle, and Fibonacci all
-- validate their inputs entirely with the generic nonnegative-integer
-- checks below. This file still exists (rather than having each
-- submodule require GABAS-LINAL's core directly) so that the moment a
-- combinatorics-specific validator IS needed, it has an obvious, already-
-- established home.
local LinalCore = require("GABAS-LINAL.modules.core")

return {
    MACHINE_EPSILON = LinalCore.MACHINE_EPSILON,
    is_finite_number = LinalCore.is_finite_number,
    is_integer = LinalCore.is_integer,
    is_positive_number = LinalCore.is_positive_number,
    is_nonneg_number = LinalCore.is_nonneg_number,
    is_positive_integer = LinalCore.is_positive_integer,
    is_nonneg_integer = LinalCore.is_nonneg_integer,
    assert_finite_number = LinalCore.assert_finite_number,
    assert_integer = LinalCore.assert_integer,
    assert_positive_number = LinalCore.assert_positive_number,
    assert_nonneg_number = LinalCore.assert_nonneg_number,
    assert_positive_integer = LinalCore.assert_positive_integer,
    assert_nonneg_integer = LinalCore.assert_nonneg_integer,
}
