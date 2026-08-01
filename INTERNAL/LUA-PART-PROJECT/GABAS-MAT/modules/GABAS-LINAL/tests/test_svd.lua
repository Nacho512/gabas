local LA = require("GABAS-LINAL.GABAS_LINAL")

local A = {{3, 0}, {0, 2}}
local U, S, V = LA.SVD(A, {verify = true})
LA.Assert_matrix_close(LA.Mat_mul(LA.Mat_mul(U, S), LA.Conjugate_transpose(V)), A, 1e-8)
local _, _, _, singular_values = LA.SVD({{1e-12, 0}, {0, 2e-12}})
LA.Assert_vector_close(singular_values, {2e-12, 1e-12}, 1e-20)
