local LA = require("GABAS-LINAL.GABAS_LINAL")

local input = {1, 2, 3, 4}
LA.Assert_vector_close(LA.IFFT(LA.FFT(input)), input, 1e-10)
LA.Assert_error(function() LA.FFT({1, 2, 3}) end, "power of 2")
