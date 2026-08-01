local LA = require("GABAS-LINAL.GABAS_LINAL")

local A = {{1, 0}, {0, 0}}
local Aplus = LA.Pinv(A, {verify = true})
LA.Assert_matrix_close(LA.Mat_mul(LA.Mat_mul(A, Aplus), A), A, 1e-8)
LA.Assert_matrix_close(LA.Mat_mul(LA.Mat_mul(Aplus, A), Aplus), Aplus, 1e-8)
