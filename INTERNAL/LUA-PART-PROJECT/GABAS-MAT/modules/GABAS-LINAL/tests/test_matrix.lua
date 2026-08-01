local LA = require("GABAS-LINAL.GABAS_LINAL")

local A = {{1, 2}, {3, 4}}
LA.Assert_matrix_close(LA.Mat_mul(A, LA.Eye(2)), A, 1e-12)
LA.Assert_vector_close(LA.Matrix_vector_mul(A, {2, 1}), {4, 10}, 1e-12)
LA.Assert_matrix_close(LA.Conjugate_transpose({{LA.Complex(1, 2)}}), {{LA.Complex(1, -2)}}, 1e-12)
LA.Assert_error(function() LA.Mat_mul({{1, 2}, {3}}, LA.Eye(2)) end, "rectangular")
