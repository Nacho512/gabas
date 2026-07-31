local LA = require("GABAS-LINAL.GABAS_LINAL")

local zero_values = LA.Eigenvalues({{0, 0}, {0, 0}})
LA.Assert_Vector_Close(zero_values, {0, 0}, 1e-12)
local values, vectors = LA.Eigenvectors({{2, 0}, {0, 3}})
for j = 1, 2 do
    local v = {vectors[1][j], vectors[2][j]}
    local Av = LA.Matrix_Vector_Mul({{2, 0}, {0, 3}}, v)
    LA.Assert_Vector_Close(Av, {values[j] * v[1], values[j] * v[2]}, 1e-8)
end
local pair = LA.Eigenvalues({{0, -1}, {1, 0}})
assert(LA.Is_Complex(pair[1]) and LA.Is_Complex(pair[2]))
LA.Assert_Error(function() LA.Eigenvectors({{0, 1}, {0, 0}}) end, "symmetric")

local systems = require("GABAS-LINAL.modules.systems")
local original_solve = systems.Solve
systems.Solve = function() error("injected solve failure") end
LA.Assert_Error(function() LA.Rayleigh({{2, 0}, {0, 3}}, {1, 1}) end, "shifted solve failed")
systems.Solve = original_solve
