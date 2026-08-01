-- [Special Functions]
-- Special functions: thin aggregator over the modules/ submodules below,
-- mirroring GABAS_CALC_ONE_VAR.lua's own aggregator pattern -- each
-- submodule's public functions get flattened into one table, so calling
-- code is unaffected by how the implementation is split across files.
--
--   special_functions -- modules/special_functions.lua -- Gamma, Log_gamma, Beta, Log_beta, Gamma_p, Gamma_q, Erf, Erfc
--
-- Claude: "core" is deliberately NOT in this list -- core.lua holds
-- internal validation plumbing that special_functions.lua depends on
-- directly, not something meant for the end user. If it were merged in
-- here too, this same loop would silently promote it to public API.
--
-- Scoped deliberately small for a first pass: Gamma and Log_gamma are
-- the foundation the rest of this family (Beta, incomplete gamma/beta,
-- erf, Bessel, ...) will build on -- each of those is its own real
-- numerical-algorithm undertaking, added incrementally, not guessed at
-- or stubbed out in bulk here.
local submodules = {
    "special_functions",
}

local SpecialFunctions = {}

for _, name in ipairs(submodules) do
    local mod = require("GABAS-SPECIAL-FUNCTIONS.modules." .. name)
    for key, value in pairs(mod) do
        assert(SpecialFunctions[key] == nil,
            "GABAS_SPECIAL_FUNCTIONS: duplicate public export '" .. key .. "' from module " .. name .. ".")
        SpecialFunctions[key] = value
    end
end

return SpecialFunctions
