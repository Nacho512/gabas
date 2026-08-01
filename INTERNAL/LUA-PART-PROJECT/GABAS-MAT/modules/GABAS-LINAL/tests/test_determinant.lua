local LA = require("GABAS-LINAL.GABAS_LINAL")

LA.Assert_close(LA.MDet({{1, 2}, {3, 4}}), -2, 1e-12)
LA.Assert_close(LA.MDet({{1, 2}, {2, 4}}), 0, 1e-12)
LA.Assert_error(function() LA.MDet({{1, 2}, {3}}) end, "rectangular")
