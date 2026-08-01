local LA = require("GABAS-LINAL.GABAS_LINAL")

local A, b = {{4, 7}, {2, 6}}, {1, 0}
LA.Assert_matrix_close(LA.Mat_mul(A, LA.Inverse(A)), LA.Eye(#A), 1e-10)
LA.Assert_vector_close(LA.Matrix_vector_mul(A, LA.Solve(A, b)), b, 1e-10)
LA.Assert_error(function() LA.Solve({{1, 2}, {3}}, b) end, "rectangular")
LA.Assert_error(function() LA.Solve(A, {1, math.huge}) end, "finite")
LA.Assert_close(LA.Inverse({{1e-12}})[1][1], 1e12, 1e-3)
