local LA = require("GABAS-LINAL.GABAS_LINAL")

local A = {{1, 0}, {0, 0}}
local Aplus = LA.Pinv(A, {verify = true})
LA.Assert_Matrix_Close(LA.Mat_Mul(LA.Mat_Mul(A, Aplus), A), A, 1e-8)
LA.Assert_Matrix_Close(LA.Mat_Mul(LA.Mat_Mul(Aplus, A), Aplus), Aplus, 1e-8)
