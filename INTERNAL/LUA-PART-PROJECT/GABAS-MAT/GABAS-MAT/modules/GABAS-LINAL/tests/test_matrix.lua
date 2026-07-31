local LA = require("GABAS-LINAL.GABAS_LINAL")

local A = {{1, 2}, {3, 4}}
LA.Assert_Matrix_Close(LA.Mat_Mul(A, LA.Eye(2)), A, 1e-12)
LA.Assert_Vector_Close(LA.Matrix_Vector_Mul(A, {2, 1}), {4, 10}, 1e-12)
LA.Assert_Matrix_Close(LA.Conjugate_Transpose({{LA.Complex(1, 2)}}), {{LA.Complex(1, -2)}}, 1e-12)
LA.Assert_Error(function() LA.Mat_Mul({{1, 2}, {3}}, LA.Eye(2)) end, "rectangular")
