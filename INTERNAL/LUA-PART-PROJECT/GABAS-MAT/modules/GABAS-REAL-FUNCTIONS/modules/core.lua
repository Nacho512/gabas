-- [Core]
-- Internal infrastructure for GABAS-REAL-FUNCTIONS's own submodules.
-- Nothing here is meant to be part of the public API -- every function
-- here is plumbing, not something an end user calls directly.
--
-- Claude: reuses GABAS-LINAL's own core.lua for every validator this
-- module needs (finite-number checks, ...) rather than re-deriving them
-- -- same pattern GABAS-CALC-ONE-VAR's, GABAS-DISCRETE-MATH's, and
-- GABAS-COMBINATORICS's own core.lua files already established.
-- Required directly from GABAS-LINAL, keeping this project an
-- independent sibling rather than needlessly coupled to any other
-- domain project.
--
-- Nothing domain-specific has been needed here yet -- every function in
-- modules/real_functions.lua validates its inputs entirely with the
-- generic finite-number check below. This file still exists (rather
-- than having that submodule require GABAS-LINAL's core directly) so
-- the moment a real-functions-specific validator IS needed, it has an
-- obvious, already-established home.
local LinalCore = require("GABAS-LINAL.modules.core")

return {
    MACHINE_EPSILON = LinalCore.MACHINE_EPSILON,
    is_finite_number = LinalCore.is_finite_number,
    assert_finite_number = LinalCore.assert_finite_number,
}
