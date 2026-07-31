local LA = require("GABAS-LINAL.GABAS_LINAL")

local A, b = {{4, 7}, {2, 6}}, {1, 0}
LA.Assert_Matrix_Close(LA.Mat_Mul(A, LA.Inverse(A)), LA.Eye(#A), 1e-10)
LA.Assert_Vector_Close(LA.Matrix_Vector_Mul(A, LA.Solve(A, b)), b, 1e-10)
LA.Assert_Error(function() LA.Solve({{1, 2}, {3}}, b) end, "rectangular")
LA.Assert_Error(function() LA.Solve(A, {1, math.huge}) end, "finite")
LA.Assert_Close(LA.Inverse({{1e-12}})[1][1], 1e12, 1e-3)
