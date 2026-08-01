-- [Real Functions]
-- Real one-variable utility functions: thin aggregator over the
-- modules/ submodules below, mirroring GABAS_CALC_ONE_VAR.lua's own
-- aggregator pattern -- each submodule's public functions get flattened
-- into one table, so calling code is unaffected by how the
-- implementation is split across files.
--
--   real_functions -- modules/real_functions.lua -- Sign, Floor, Ceil, Round, Trunc, Frac, Min, Max, Clamp, Hypot
--
-- Claude: "core" is deliberately NOT in this list -- core.lua holds
-- internal validation plumbing that real_functions.lua depends on
-- directly, not something meant for the end user. If it were merged in
-- here too, this same loop would silently promote it to public API.
--
-- Sin/cos/tan/exp/ln/sqrt/abs and their inverse/hyperbolic relatives are
-- NOT here (yet) -- see real_functions.lua's own header for exactly
-- why (they already live, properly, as Dual/Complex-aware dispatchers
-- in GABAS-CALC-ONE-VAR/modules/compile_expression.lua, and moving them
-- here is a deliberate future step, not this one).
local submodules = {
    "real_functions",
}

local RealFunctions = {}

for _, name in ipairs(submodules) do
    local mod = require("GABAS-REAL-FUNCTIONS.modules." .. name)
    for key, value in pairs(mod) do
        assert(RealFunctions[key] == nil,
            "GABAS_REAL_FUNCTIONS: duplicate public export '" .. key .. "' from module " .. name .. ".")
        RealFunctions[key] = value
    end
end

return RealFunctions
