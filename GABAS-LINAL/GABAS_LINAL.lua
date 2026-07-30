-- ===========================================================================
-- GABAS_LINAL: thin aggregator over the modules/ submodules below. Splitting
-- the actual implementation into per-category files keeps each one small
-- enough to read and edit in isolation; this file's only job is to require
-- every domain submodule and flatten their public functions into one table,
-- so calling code is completely unaffected: `GL.MDet(A)`, `GL.FFT(v)`, etc.
-- work exactly as they did when everything lived in a single file.
--
--   Core        -- modules/core.lua        -- internal plumbing, not public
--   Complex     -- modules/complex.lua
--   Matrix      -- modules/matrix.lua
--   Tensor      -- modules/tensor.lua
--   Vector      -- modules/vector.lua
--   Determinant -- modules/determinant.lua
--   Systems     -- modules/systems.lua
--   Eigen       -- modules/eigen.lua
--   FFT         -- modules/fft.lua
-- ===========================================================================

local GL = {}

-- "core" is deliberately NOT in this list: modules/core.lua holds internal
-- plumbing (EPSILON, cabs, copy_matrix, find_pivot_row, next_pow2) that the
-- domain submodules below require directly from each other, not something
-- meant for the end user. If it were merged in here too, this same loop
-- would silently promote it to public API.
local submodules = {
    "complex",
    "matrix",
    "tensor",
    "vector",
    "determinant",
    "systems",
    "eigen",
    "fft",
}

for _, name in ipairs(submodules) do
    local mod = require("GABAS-LINAL.modules." .. name)
    for key, value in pairs(mod) do
        GL[key] = value
    end
end

return GL
