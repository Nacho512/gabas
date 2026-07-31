local LA = require("GABAS-LINAL.GABAS_LINAL")

local tensor = LA.Tensor({1, 2, 3, 4, 5, 6, 7, 8}, {2, 2, 2})
assert(tensor[2][2][2] == 8)
local eye = LA.Eye_Tensor({2, 2, 2})
assert(eye[1][1][1] == 1 and eye[1][1][2] == 0 and eye[2][2][2] == 1)
LA.Assert_Error(function() LA.Tensor({1, math.huge}, {1, 1, 2}) end, "finite")
