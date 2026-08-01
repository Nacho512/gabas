-- [Calc One Var]
-- One-variable calculus: thin aggregator over the modules/ submodules
-- below, mirroring GABAS_LINAL.lua's own aggregator pattern -- each
-- submodule's public functions get flattened into one table, so calling
-- code is unaffected by how the implementation is split across files.
--
--   compile_expression   -- modules/compile_expression.lua -- parses/compiles a symbolic expression into a callable Lua function
--   dif_numeric          -- modules/dif_numeric.lua         -- numerical differentiation (forward-mode AD via Dual numbers)
--   integration_numeric  -- modules/integration_numeric.lua -- numerical integration (adaptive Gauss-Kronrod quadrature)
--   roots_numeric        -- modules/roots_numeric.lua       -- numerical root-finding (Bisection, Newton, Secant)
--
-- Claude: "core", "dual", and "quadrature" are deliberately NOT in this
-- list -- core.lua holds internal validation plumbing, dual.lua and
-- quadrature.lua hold internal numerical building blocks that
-- dif_numeric.lua/integration_numeric.lua depend on directly. None of the
-- three are meant to be part of the public API; if they were merged in
-- here too, this same loop would silently promote them to it.
local submodules = {
    "compile_expression",
    "dif_numeric",
    "integration_numeric",
    "roots_numeric",
}

local CalcOneVar = {}

for _, name in ipairs(submodules) do
    local mod = require("GABAS-CALC-ONE-VAR.modules." .. name)
    for key, value in pairs(mod) do
        assert(CalcOneVar[key] == nil,
            "GABAS_CALC_ONE_VAR: duplicate public export '" .. key .. "' from module " .. name .. ".")
        CalcOneVar[key] = value
    end
end

return CalcOneVar
