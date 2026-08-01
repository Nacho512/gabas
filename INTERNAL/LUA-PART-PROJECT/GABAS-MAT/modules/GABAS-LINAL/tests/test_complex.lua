local LA = require("GABAS-LINAL.GABAS_LINAL")

local z = LA.Complex(3, 4)
LA.Assert_close(z:abs(), 5, 1e-12)
LA.Assert_close((LA.Complex(1, 0) / LA.Complex(1e-12, 0)).re, 1e12, 1e-3)
LA.Assert_close((LA.Complex(1, 0) / LA.Complex(1e-300, 0)).re, 1e300, 1e286)
LA.Assert_error(function() return z / LA.Complex(0, 0) end, "division by zero")
LA.Assert_error(function() return LA.Complex(math.huge, 0) end, "finite numbers")
assert(LA.Is_complex(z))

-- Asin/Acos/Atan: real values, complex round-trip, and the 0^v edge case
-- (Sqrt(0) = 0^0.5) that Asin(1) surfaces via its own formula.
LA.Assert_close(LA.Asin(1).re, math.pi / 2, 1e-12)
LA.Assert_close(LA.Acos(1).re, 0, 1e-12)
LA.Assert_close(LA.Atan(1).re, math.pi / 4, 1e-12)
local w = LA.Complex(0.3, 0.7)
LA.Assert_vector_close({LA.Sin(LA.Asin(w)).re, LA.Sin(LA.Asin(w)).im}, {w.re, w.im}, 1e-12)
LA.Assert_vector_close({LA.Cos(LA.Acos(w)).re, LA.Cos(LA.Acos(w)).im}, {w.re, w.im}, 1e-12)
LA.Assert_vector_close({LA.Tan(LA.Atan(w)).re, LA.Tan(LA.Atan(w)).im}, {w.re, w.im}, 1e-12)
LA.Assert_close((LA.Complex(0, 0) ^ 0.5).re, 0, 1e-12)
LA.Assert_error(function() return LA.Complex(0, 0) ^ (-0.5) end, "nonpositive real part")
