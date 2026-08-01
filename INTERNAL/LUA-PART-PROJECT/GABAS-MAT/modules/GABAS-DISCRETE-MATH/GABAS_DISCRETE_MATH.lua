-- [Discrete Math]
-- Discrete mathematics: thin aggregator over the modules/ submodules
-- below, mirroring GABAS_CALC_ONE_VAR.lua's own aggregator pattern --
-- each submodule's public functions get flattened into one table, so
-- calling code is unaffected by how the implementation is split across
-- files.
--
--   discrete_math -- modules/discrete_math.lua -- numeral-base conversion, GCD/LCM, primality
--
-- Claude: "core" is deliberately NOT in this list -- core.lua holds
-- internal validation plumbing that discrete_math.lua depends on
-- directly, not something meant for the end user. If it were merged in
-- here too, this same loop would silently promote it to public API.
local submodules = {
    "discrete_math",
}

local DiscreteMath = {}

for _, name in ipairs(submodules) do
    local mod = require("GABAS-DISCRETE-MATH.modules." .. name)
    for key, value in pairs(mod) do
        assert(DiscreteMath[key] == nil,
            "GABAS_DISCRETE_MATH: duplicate public export '" .. key .. "' from module " .. name .. ".")
        DiscreteMath[key] = value
    end
end

return DiscreteMath
