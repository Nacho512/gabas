local CalcOneVar = require("GABAS-CALC-ONE-VAR.GABAS_CALC_ONE_VAR")
local Testing = require("GABAS-LINAL.modules.testing")

local Inf = CalcOneVar.Integral_infinite

-- Semi-infinite intervals (technical document, section 34).
Testing.Assert_close((Inf(function(x) return math.exp(-x) end, 0, math.huge)), 1, 1e-6,
    "integral of e^-x on [0,inf)")
Testing.Assert_close((Inf(function(x) return math.exp(x) end, -math.huge, 0)), 1, 1e-6,
    "integral of e^x on (-inf,0]")
Testing.Assert_close((Inf(function(x) return 1 / (1 + x * x) end, 0, math.huge)), math.pi / 2, 1e-6,
    "integral of 1/(1+x^2) on [0,inf)")

-- Whole real line.
Testing.Assert_close((Inf(function(x) return math.exp(-x * x) end, -math.huge, math.huge)), math.sqrt(math.pi), 1e-6,
    "integral of e^(-x^2) on (-inf,inf)")
Testing.Assert_close((Inf(function(x) return 1 / (1 + x * x) end, -math.huge, math.huge)), math.pi, 1e-6,
    "integral of 1/(1+x^2) on (-inf,inf)")

-- Both limits already finite: delegates straight to Integral_definite,
-- no transform needed.
Testing.Assert_close((Inf(math.sin, 0, math.pi)), 2, 1e-8, "delegates when both limits are finite")

-- Fast-decaying power-law tail (p comfortably above the convergence
-- boundary) still converges cleanly.
local v19, e19, s19 = Inf(function(x) return x ^ (-1.9) end, 1, math.huge, {max_intervals = 300})
Testing.Assert_close(v19, 1 / 0.9, 1e-4, "x^-1.9 on [1,inf)")
assert(s19 == "SUCCESS", "x^-1.9 should converge cleanly")

-- Known, documented limitation: as p approaches the convergence boundary
-- (p=1) from above, the algebraic transform turns the slowly-decaying
-- tail into a genuine endpoint singularity in the transformed variable
-- (verified: x^-p on [1,inf) transforms to (1-t)^(p-2) on [0,1], which
-- is singular at t=1 whenever p<2). This ordinary adaptive rule has no
-- singularity-specific handling, so it must not report a false SUCCESS
-- once p gets close enough to 1 -- tanh-sinh (not yet built) is the
-- documented fallback for this combination.
local _, _, s15 = Inf(function(x) return x ^ (-1.5) end, 1, math.huge, {max_intervals = 300})
assert(s15 ~= "SUCCESS", "x^-1.5 on [1,inf) must not report false convergence")

-- Divergent integral (p=1, harmonic-like tail): must not report SUCCESS.
local _, _, s1 = Inf(function(x) return 1 / x end, 1, math.huge, {max_intervals = 200})
assert(s1 ~= "SUCCESS", "the divergent integral of 1/x on [1,inf) must not report SUCCESS")

-- Input validation.
Testing.Assert_error(function() return Inf(math.sin, math.huge, math.huge) end, "must be less than")
Testing.Assert_error(function() return Inf(42, 0, math.huge) end, "Integral_infinite: f")
