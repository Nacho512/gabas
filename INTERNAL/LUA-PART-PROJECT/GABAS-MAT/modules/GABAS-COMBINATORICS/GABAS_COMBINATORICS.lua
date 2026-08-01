-- [Combinatorics]
-- Combinatorics: thin aggregator over the modules/ submodules below,
-- mirroring GABAS_CALC_ONE_VAR.lua's own aggregator pattern -- each
-- submodule's public functions get flattened into one table, so calling
-- code is unaffected by how the implementation is split across files.
--
--   combinatorics -- modules/combinatorics.lua -- Factorial, Binomial_coefficient, Pascal's triangle, Fibonacci
--
-- Claude: "core" is deliberately NOT in this list -- core.lua holds
-- internal validation plumbing that combinatorics.lua depends on
-- directly, not something meant for the end user. If it were merged in
-- here too, this same loop would silently promote it to public API.
local submodules = {
    "combinatorics",
}

local Combinatorics = {}

for _, name in ipairs(submodules) do
    local mod = require("GABAS-COMBINATORICS.modules." .. name)
    for key, value in pairs(mod) do
        assert(Combinatorics[key] == nil,
            "GABAS_COMBINATORICS: duplicate public export '" .. key .. "' from module " .. name .. ".")
        Combinatorics[key] = value
    end
end

return Combinatorics
