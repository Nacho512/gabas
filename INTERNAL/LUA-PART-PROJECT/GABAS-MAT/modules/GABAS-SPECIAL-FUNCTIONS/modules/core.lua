-- [Core]
-- Internal infrastructure for GABAS-SPECIAL-FUNCTIONS's own submodules.
-- Nothing here is meant to be part of the public API -- every function
-- here is plumbing, not something an end user calls directly.
--
-- Claude: reuses GABAS-LINAL's own core.lua for every validator this
-- module needs (finite-number checks, ...) rather than re-deriving them
-- -- same pattern GABAS-CALC-ONE-VAR's, GABAS-DISCRETE-MATH's,
-- GABAS-COMBINATORICS's, and GABAS-REAL-FUNCTIONS's own core.lua files
-- already established. Required directly from GABAS-LINAL, keeping this
-- project an independent sibling.
local LinalCore = require("GABAS-LINAL.modules.core")

return {
    MACHINE_EPSILON = LinalCore.MACHINE_EPSILON,
    is_finite_number = LinalCore.is_finite_number,
    assert_finite_number = LinalCore.assert_finite_number,
}
