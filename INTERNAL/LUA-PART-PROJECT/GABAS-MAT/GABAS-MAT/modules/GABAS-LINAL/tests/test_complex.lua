local LA = require("GABAS-LINAL.GABAS_LINAL")

local z = LA.Complex(3, 4)
LA.Assert_Close(z:abs(), 5, 1e-12)
LA.Assert_Close((LA.Complex(1, 0) / LA.Complex(1e-12, 0)).re, 1e12, 1e-3)
LA.Assert_Close((LA.Complex(1, 0) / LA.Complex(1e-300, 0)).re, 1e300, 1e286)
LA.Assert_Error(function() return z / LA.Complex(0, 0) end, "division by zero")
LA.Assert_Error(function() return LA.Complex(math.huge, 0) end, "finite numbers")
assert(LA.Is_Complex(z))
