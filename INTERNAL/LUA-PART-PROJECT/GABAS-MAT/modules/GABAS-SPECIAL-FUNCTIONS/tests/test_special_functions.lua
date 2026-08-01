local SpecialFunctions = require("GABAS-SPECIAL-FUNCTIONS.GABAS_SPECIAL_FUNCTIONS")
local Combinatorics = require("GABAS-COMBINATORICS.GABAS_COMBINATORICS")
local Testing = require("GABAS-LINAL.modules.testing")

-- ===== Gamma: known closed-form and integer values =====

-- Gamma(n+1) = n! for every nonnegative integer n -- cross-checked
-- directly against GABAS-COMBINATORICS's own, independently verified
-- Factorial (differential testing, not just internal self-consistency).
for n = 0, 10 do
    Testing.Assert_close(SpecialFunctions.Gamma(n + 1), Combinatorics.Factorial(n), 1e-6,
        "Gamma(" .. (n + 1) .. ") == " .. n .. "!")
end

-- Gamma(1/2) = sqrt(pi) (the classical closed form).
Testing.Assert_close(SpecialFunctions.Gamma(0.5), math.sqrt(math.pi), 1e-12, "Gamma(0.5) = sqrt(pi)")

-- Gamma(-1/2) = -2*sqrt(pi) -- exercises the reflection formula AND its
-- sign handling for negative x.
Testing.Assert_close(SpecialFunctions.Gamma(-0.5), -2 * math.sqrt(math.pi), 1e-12, "Gamma(-0.5) = -2*sqrt(pi)")

-- Reference values below (mpmath-independent: Python's own math.gamma)
-- for non-integer, non-half-integer x.
Testing.Assert_close(SpecialFunctions.Gamma(0.3), 2.991568987687591, 1e-9, "Gamma(0.3)")
Testing.Assert_close(SpecialFunctions.Gamma(3.7), 4.170651783796603, 1e-9, "Gamma(3.7)")
Testing.Assert_close(SpecialFunctions.Gamma(-4.5), -0.06001960130050425, 1e-9, "Gamma(-4.5)")

-- Poles at every nonpositive integer -- must be rejected explicitly, not
-- silently returned as Inf/NaN.
Testing.Assert_error(function() return SpecialFunctions.Gamma(0) end, "Gamma: x must not be a nonpositive integer")
Testing.Assert_error(function() return SpecialFunctions.Gamma(-1) end, "Gamma: x must not be a nonpositive integer")
Testing.Assert_error(function() return SpecialFunctions.Gamma(-7) end, "Gamma: x must not be a nonpositive integer")
Testing.Assert_error(function() return SpecialFunctions.Gamma("2") end, "Gamma: x")

-- ===== Log_gamma =====
--
-- Reference values below are Python's math.lgamma (an independent,
-- trusted reference), not values this implementation itself produced.

Testing.Assert_close(SpecialFunctions.Log_gamma(100), 359.1342053695754, 1e-9, "Log_gamma(100)")
Testing.Assert_close(SpecialFunctions.Log_gamma(170), 701.437263808737, 1e-9, "Log_gamma(170)")
Testing.Assert_close(SpecialFunctions.Log_gamma(0.3), 1.0957979948180752, 1e-9, "Log_gamma(0.3)")
Testing.Assert_close(SpecialFunctions.Log_gamma(-4.5), -2.8130840817693166, 1e-9, "Log_gamma(-4.5)")

-- Internal consistency: Log_gamma(x) == log(|Gamma(x)|) for moderate x
-- (where Gamma(x) itself doesn't overflow).
for _, x in ipairs({1.2, 2.5, 4.4, -0.5, -3.3, 7}) do
    Testing.Assert_close(SpecialFunctions.Log_gamma(x), math.log(math.abs(SpecialFunctions.Gamma(x))), 1e-8,
        "Log_gamma(x) == log(|Gamma(x)|) at x=" .. x)
end

Testing.Assert_error(function() return SpecialFunctions.Log_gamma(-2) end,
    "Log_gamma: x must not be a nonpositive integer")

-- ===== The overflow boundary: Gamma legitimately overflows to math.huge
-- once the true value exceeds a double's range (~Gamma(172) and up),
-- while Log_gamma -- staying in log-space the whole time -- does not. =====

assert(SpecialFunctions.Gamma(200) == math.huge, "Gamma(200) must overflow to +inf (the true value exceeds a double)")
Testing.Assert_close(SpecialFunctions.Log_gamma(200), 857.9336698258575, 1e-9,
    "Log_gamma(200) stays finite and accurate even where Gamma itself overflows")

-- ===== Beta / Log_beta =====
--
-- Reference values below are mpmath's own beta() at 30-digit working
-- precision (an independent, trusted reference distinct from this
-- module's own Log_gamma-based formula), plus a cross-check against
-- Beta's DEFINING INTEGRAL for a few values -- two independent
-- verifications, not just one formula checked against itself.

Testing.Assert_close(SpecialFunctions.Beta(1, 1), 1, 1e-12, "Beta(1,1) = 1")
Testing.Assert_close(SpecialFunctions.Beta(0.5, 0.5), math.pi, 1e-9, "Beta(0.5,0.5) = pi")
Testing.Assert_close(SpecialFunctions.Beta(2, 3), 1 / 12, 1e-12, "Beta(2,3) = 1/12")
Testing.Assert_close(SpecialFunctions.Beta(10, 20), 4.9925087406346777e-09, 1e-9, "Beta(10,20)")
Testing.Assert_close(SpecialFunctions.Beta(2.5, 3.7), 0.032727368606257838, 1e-9, "Beta(2.5,3.7)")
Testing.Assert_close(SpecialFunctions.Beta(0.1, 0.1), 19.714639489050161, 1e-8, "Beta(0.1,0.1)")

-- The overflow case: Gamma(100)*Gamma(100)/Gamma(200) computed naively
-- would divide by +inf (Gamma(200) already overflows a double) and give
-- a silently WRONG 0 -- Beta(100,100) is actually an ordinary, tiny but
-- nonzero positive number. Verified this really is what the naive
-- formula does, in Python, before writing Log_beta this way.
Testing.Assert_close(SpecialFunctions.Beta(100, 100), 2.2087606931995026e-61, 1e-70, "Beta(100,100), the overflow case")

-- Cross-check against the defining integral, via a fine Riemann sum
-- (independent of Log_gamma entirely) -- not exact quadrature, so a
-- looser tolerance.
local function beta_by_integral(a, b, n)
    local sum, h = 0, 1 / n
    for i = 1, n - 1 do
        local t = i * h
        sum = sum + (t ^ (a - 1)) * ((1 - t) ^ (b - 1))
    end
    return sum * h
end
Testing.Assert_close(SpecialFunctions.Beta(2, 3), beta_by_integral(2, 3, 200000), 1e-4,
    "Beta(2,3) matches its own defining integral, computed independently")
Testing.Assert_close(SpecialFunctions.Beta(5, 5), beta_by_integral(5, 5, 200000), 1e-6,
    "Beta(5,5) matches its own defining integral, computed independently")

-- Log_beta internal consistency.
for _, ab in ipairs({{1, 1}, {2, 3}, {10, 20}, {100, 100}}) do
    Testing.Assert_close(SpecialFunctions.Log_beta(ab[1], ab[2]), math.log(SpecialFunctions.Beta(ab[1], ab[2])), 1e-8,
        "Log_beta == log(Beta) at a=" .. ab[1] .. ",b=" .. ab[2])
end

-- Input validation: a, b must both be strictly positive.
Testing.Assert_error(function() return SpecialFunctions.Beta(0, 1) end, "Beta: a must be positive")
Testing.Assert_error(function() return SpecialFunctions.Beta(1, -1) end, "Beta: b must be positive")
Testing.Assert_error(function() return SpecialFunctions.Beta("1", 1) end, "Beta: a")
Testing.Assert_error(function() return SpecialFunctions.Log_beta(0, 1) end, "Log_beta: a must be positive")

-- ===== Gamma_p / Gamma_q =====

Testing.Assert_close(SpecialFunctions.Gamma_p(1, 1), 1 - math.exp(-1), 1e-12, "Gamma_p(1,x) = 1-e^-x (exponential CDF)")
assert(SpecialFunctions.Gamma_p(2, 0) == 0, "Gamma_p(a,0) == 0")
assert(SpecialFunctions.Gamma_q(2, 0) == 1, "Gamma_q(a,0) == 1")
for _, ax in ipairs({{1, 1}, {2, 5}, {0.5, 3}, {10, 8}, {10, 15}}) do
    local p, q = SpecialFunctions.Gamma_p(ax[1], ax[2]), SpecialFunctions.Gamma_q(ax[1], ax[2])
    Testing.Assert_close(p + q, 1, 1e-9, "Gamma_p+Gamma_q == 1 at a=" .. ax[1] .. ",x=" .. ax[2])
end

Testing.Assert_error(function() return SpecialFunctions.Gamma_p(0, 1) end, "Gamma_p: a must be positive")
Testing.Assert_error(function() return SpecialFunctions.Gamma_p(1, -1) end, "Gamma_p: x must be nonnegative")

-- ===== Erf / Erfc =====
--
-- Reference values below are Python's math.erf/math.erfc (independent,
-- trusted references), not values this implementation itself produced.

local erf_ref = {
    [0.001] = 0.0011283787909692365, [0.1] = 0.1124629160182849, [0.5] = 0.5204998778130465,
    [1] = 0.8427007929497149, [1.5] = 0.9661051464753108, [2] = 0.9953222650189527,
    [2.5] = 0.999593047982555, [3] = 0.9999779095030014, [4] = 0.9999999845827421,
    [-0.5] = -0.5204998778130465, [-1] = -0.8427007929497149, [-2] = -0.9953222650189527,
}
for x, ref in pairs(erf_ref) do
    Testing.Assert_close(SpecialFunctions.Erf(x), ref, 1e-9, "Erf(" .. x .. ")")
end
assert(SpecialFunctions.Erf(0) == 0, "Erf(0) == 0")

local erfc_ref = {
    [0.001] = 0.9988716212090307, [0.1] = 0.8875370839817152, [0.5] = 0.4795001221869535,
    [1] = 0.15729920705028513, [2] = 0.004677734981047265, [3] = 2.2090496998585438e-05,
    [4] = 1.541725790028002e-08, [5] = 1.5374597944280351e-12, [6] = 2.1519736712498916e-17,
    [-0.5] = 1.5204998778130465, [-1] = 1.842700792949715, [-2] = 1.9953222650189528,
}
for x, ref in pairs(erfc_ref) do
    Testing.Assert_close(SpecialFunctions.Erfc(x), ref, math.max(1e-9, math.abs(ref) * 1e-6), "Erfc(" .. x .. ")")
end

-- The catastrophic-cancellation case this whole design avoids: for large
-- x, Erf(x) itself rounds to exactly 1.0, so "1 - Erf(x)" would silently
-- give exactly 0 -- but the true Erfc(x) keeps shrinking meaningfully
-- for a long time after. Verified directly against Python out to x=20.
assert(SpecialFunctions.Erf(10) == 1, "Erf(10) rounds to exactly 1.0 in double precision")
Testing.Assert_close(SpecialFunctions.Erfc(10), 2.088487583762545e-45, 1e-55,
    "Erfc(10) stays accurate even though 1-Erf(10) would silently give 0")
Testing.Assert_close(SpecialFunctions.Erfc(15), 7.212994172451206e-100, 1e-110,
    "Erfc(15) stays accurate")
Testing.Assert_close(SpecialFunctions.Erfc(20), 5.3958656116079005e-176, 1e-186,
    "Erfc(20) stays accurate")

Testing.Assert_error(function() return SpecialFunctions.Erf("1") end, "Erf: x")
Testing.Assert_error(function() return SpecialFunctions.Erfc("1") end, "Erfc: x")

-- ===== Bessel_j =====
--
-- Reference values below are mpmath's besselj() at 30-digit working
-- precision (an independent, trusted reference), not values this
-- implementation itself produced.

local bessel_j_ref = {
    {0, 0.5, 0.9384698072408129}, {0, 1, 0.76519768655796655}, {0, 5, -0.1775967713143383},
    {0, 10, -0.24593576445134834}, {0, 15, -0.014224472826780773},
    {1, 0.5, 0.24226845767487389}, {1, 1, 0.44005058574493352}, {1, 5, -0.32757913759146522},
    {1, 10, 0.043472746168861437}, {1, 15, 0.20510403861352276},
    {2, 3, 0.48609126058589108}, {2, 11, 0.13904751877870127},
    {3, 7, -0.16755558799533424}, {5, 13, 0.13161955992748079},
}
for _, case in ipairs(bessel_j_ref) do
    local n, x, ref = case[1], case[2], case[3]
    Testing.Assert_close(SpecialFunctions.Bessel_j(n, x), ref, 1e-9, "Bessel_j(" .. n .. "," .. x .. ")")
end

-- Beyond the old power series' |x| <= 15 boundary: Miller's algorithm
-- (backward recurrence + Neumann normalization) has no known breakdown
-- mechanism as |x| or n grow, and this is verified directly against
-- mpmath.besselj across large x, large n, n >> x, and x >> n together --
-- not just a couple of spot checks past the old boundary.
local bessel_j_extended_ref = {
    {0, 20, 0.16702466434058315}, {5, 20, 0.15116976798239497}, {20, 20, 0.16474777377532653},
    {0, 50, 0.055812327669251815}, {10, 50, -0.11384784914946939}, {50, 50, 0.12140902189761506},
    {0, 100, 0.019985850304223122}, {25, 100, 0.078504273355993287}, {100, 100, 0.096366673295861560},
    {50, 5, 2.2942476159525401e-45}, {100, 1, 8.4318287896267085e-189},
    {200, 300, -0.019369872600834379},
}
for _, case in ipairs(bessel_j_extended_ref) do
    local n, x, ref = case[1], case[2], case[3]
    local tol = math.max(math.abs(ref) * 1e-10, 1e-55)
    Testing.Assert_close(SpecialFunctions.Bessel_j(n, x), ref, tol, "Bessel_j(" .. n .. "," .. x .. ") [extended domain]")
end

-- Negative x beyond the old boundary too: J_n(-x) = (-1)^n * J_n(x) for
-- integer n, so this also exercises the extended-domain path both
-- directly (mpmath.besselj at negative x) and via the identity.
local bessel_j_negative_extended_ref = {
    {0, -20, 0.16702466434058315}, {1, -20, -0.066833124175850046}, {3, -50, -0.092734804061634432},
}
for _, case in ipairs(bessel_j_negative_extended_ref) do
    local n, x, ref = case[1], case[2], case[3]
    Testing.Assert_close(SpecialFunctions.Bessel_j(n, x), ref, math.abs(ref) * 1e-10, "Bessel_j(" .. n .. "," .. x .. ") [extended, negative x]")
end

-- J_n(0): 1 for n=0, 0 for n>0.
assert(SpecialFunctions.Bessel_j(0, 0) == 1, "Bessel_j(0,0) == 1")
assert(SpecialFunctions.Bessel_j(3, 0) == 0, "Bessel_j(n>0,0) == 0")

-- J_n(-x) = (-1)^n * J_n(x).
for _, n in ipairs({0, 1, 2, 3}) do
    local x = 7.3
    local expected = ((n % 2 == 0) and 1 or -1) * SpecialFunctions.Bessel_j(n, x)
    Testing.Assert_close(SpecialFunctions.Bessel_j(n, -x), expected, 1e-9,
        "Bessel_j(" .. n .. ",-x) = (-1)^" .. n .. "*Bessel_j(" .. n .. ",x)")
end

-- Input validation.
Testing.Assert_error(function() return SpecialFunctions.Bessel_j(-1, 1) end, "Bessel_j: n must be a nonnegative integer")
Testing.Assert_error(function() return SpecialFunctions.Bessel_j(1.5, 1) end, "Bessel_j: n must be a nonnegative integer")
Testing.Assert_error(function() return SpecialFunctions.Bessel_j(0, "1") end, "Bessel_j: x")

-- ===== Bessel_y =====
--
-- Reference values below are mpmath's bessely() at 30-digit working
-- precision (an independent, trusted reference), not values this
-- implementation itself produced.

local bessel_y_ref = {
    {0, 0.1, -1.5342386513503668}, {0, 0.5, -0.44451873350670656}, {0, 1, 0.088256964215676958},
    {0, 5, -0.30851762524903378}, {0, 10, 0.055671167283599391}, {0, 12, -0.22523731263436143},
    {1, 0.1, -6.4589510947020266}, {1, 0.5, -1.4714723926702431}, {1, 1, -0.78121282130028872},
    {1, 5, 0.14786314339122684}, {1, 10, 0.24901542420695388}, {1, 12, -0.057099218260896521},
    {2, 1, -1.6506826068162544}, {2, 5, 0.36766288260552452}, {2, 10, -0.0058680824422086146},
    {3, 5, 0.14626716269319277}, {5, 8, 0.25640106499011348}, {10, 10, -0.35981415218340272},
    {20, 12, -79.349697401970764},
}
for _, case in ipairs(bessel_y_ref) do
    local n, x, ref = case[1], case[2], case[3]
    Testing.Assert_close(SpecialFunctions.Bessel_y(n, x), ref, math.abs(ref) * 1e-9, "Bessel_y(" .. n .. "," .. x .. ")")
end

-- Input validation: x=0 (logarithmic singularity) and x<0 (not
-- real-valued for Y_n) are both rejected explicitly, as is the verified
-- x <= 12 domain boundary -- the same "enforce the checked domain, don't
-- silently extrapolate" discipline Bessel_j used before Miller's
-- algorithm removed the need for it there.
Testing.Assert_error(function() return SpecialFunctions.Bessel_y(0, 0) end, "Bessel_y: x must be positive")
Testing.Assert_error(function() return SpecialFunctions.Bessel_y(0, -1) end, "Bessel_y: x must be positive")
Testing.Assert_error(function() return SpecialFunctions.Bessel_y(0, 13) end, "Bessel_y: x must be <= 12")
Testing.Assert_error(function() return SpecialFunctions.Bessel_y(-1, 1) end, "Bessel_y: n must be a nonnegative integer")
Testing.Assert_error(function() return SpecialFunctions.Bessel_y(1.5, 1) end, "Bessel_y: n must be a nonnegative integer")
Testing.Assert_error(function() return SpecialFunctions.Bessel_y(0, "1") end, "Bessel_y: x")
