local SpecialFunctions = require("GABAS-SPECIAL-FUNCTIONS.GABAS_SPECIAL_FUNCTIONS")
local Testing = require("GABAS-LINAL.modules.testing")

-- ===== Carlson_rc =====
--
-- Reference values below are mpmath's elliprc() at 30-digit working
-- precision (an independent, trusted reference), not values this
-- implementation itself produced.

local carlson_rc_ref = {
    {5, 10, 0.35124073655203632}, {0, 5, 0.70248147310407264}, {2, 2, 0.70710678118654752},
    {0.001, 3, 0.89650760127549780}, {100, 0.01, 0.52985573000990173},
}
for _, case in ipairs(carlson_rc_ref) do
    local x, y, ref = case[1], case[2], case[3]
    Testing.Assert_close(SpecialFunctions.Carlson_rc(x, y), ref, math.abs(ref) * 1e-12, "Carlson_rc(" .. x .. "," .. y .. ")")
end

-- Identity (the reference document's own test M7): RC(a,a) = 1/sqrt(a).
for _, a in ipairs({0.5, 1, 3.7, 100}) do
    Testing.Assert_close(SpecialFunctions.Carlson_rc(a, a), 1 / math.sqrt(a), 1e-12, "Carlson_rc(" .. a .. "," .. a .. ") = 1/sqrt(a)")
end

Testing.Assert_error(function() return SpecialFunctions.Carlson_rc(-1, 5) end, "Carlson_rc: x must be nonnegative")
Testing.Assert_error(function() return SpecialFunctions.Carlson_rc(1, 0) end, "Carlson_rc: y must be positive")
Testing.Assert_error(function() return SpecialFunctions.Carlson_rc(1, -1) end, "Carlson_rc: y must be positive")
Testing.Assert_error(function() return SpecialFunctions.Carlson_rc("1", 1) end, "Carlson_rc: x")

-- ===== Carlson_rf =====
--
-- Reference values below are mpmath's elliprf() at 30-digit working
-- precision.

local carlson_rf_ref = {
    {1, 2, 3, 0.72694593546890820}, {0, 2, 3, 1.0010773804561062}, {0, 1, 1, 1.5707963267948966},
    {5, 5, 5, 0.44721359549995794}, {0.001, 1, 10, 0.80543540693082469},
}
for _, case in ipairs(carlson_rf_ref) do
    local x, y, z, ref = case[1], case[2], case[3], case[4]
    Testing.Assert_close(SpecialFunctions.Carlson_rf(x, y, z), ref, math.abs(ref) * 1e-12, "Carlson_rf(" .. x .. "," .. y .. "," .. z .. ")")
end

-- Identities from the reference document's own validation plan:
-- M4 (RF invariant under permutation) and equal-argument closed form
-- RF(a,a,a) = 1/sqrt(a).
for _, a in ipairs({0.5, 1, 3.7, 100}) do
    Testing.Assert_close(SpecialFunctions.Carlson_rf(a, a, a), 1 / math.sqrt(a), 1e-12, "Carlson_rf(" .. a .. "," .. a .. "," .. a .. ") = 1/sqrt(a)")
end
do
    local x, y, z = 2, 5, 11
    local base = SpecialFunctions.Carlson_rf(x, y, z)
    local perms = {{x, z, y}, {y, x, z}, {y, z, x}, {z, x, y}, {z, y, x}}
    for _, perm in ipairs(perms) do
        Testing.Assert_close(SpecialFunctions.Carlson_rf(perm[1], perm[2], perm[3]), base, 1e-12,
            "Carlson_rf is permutation-invariant")
    end
end

-- M6: homogeneity, RF(lambda*x,lambda*y,lambda*z) = RF(x,y,z)/sqrt(lambda).
do
    local x, y, z = 3, 7, 15
    local base = SpecialFunctions.Carlson_rf(x, y, z)
    for _, lambda in ipairs({4, 1024, 1 / 16}) do
        Testing.Assert_close(SpecialFunctions.Carlson_rf(lambda * x, lambda * y, lambda * z), base / math.sqrt(lambda), 1e-11,
            "Carlson_rf homogeneity at lambda=" .. lambda)
    end
end

-- M7: RC(x,y) = RF(x,y,y).
for _, case in ipairs({{5, 10}, {0.5, 3}, {20, 1}}) do
    local x, y = case[1], case[2]
    Testing.Assert_close(SpecialFunctions.Carlson_rc(x, y), SpecialFunctions.Carlson_rf(x, y, y), 1e-12,
        "Carlson_rc(x,y) == Carlson_rf(x,y,y)")
end

Testing.Assert_error(function() return SpecialFunctions.Carlson_rf(-1, 1, 1) end, "Carlson_rf: x, y, and z must all be nonnegative")
Testing.Assert_error(function() return SpecialFunctions.Carlson_rf(0, 0, 1) end, "Carlson_rf: at most one of x, y, z may be zero")
Testing.Assert_error(function() return SpecialFunctions.Carlson_rf("1", 1, 1) end, "Carlson_rf: x")

-- ===== Carlson_rj =====
--
-- Reference values below are mpmath's elliprj() at 30-digit working
-- precision.

local carlson_rj_ref = {
    {1, 2, 3, 4, 0.23984809974956776}, {0, 1, 2, 3, 0.77688623778582332}, {2, 2, 2, 2, 0.35355339059327376},
    {0.001, 1, 10, 5, 0.26551788124327310}, {1, 1, 1, 1, 1.0},
}
for _, case in ipairs(carlson_rj_ref) do
    local x, y, z, p, ref = case[1], case[2], case[3], case[4], case[5]
    local tol = math.max(math.abs(ref) * 1e-11, 1e-12)
    Testing.Assert_close(SpecialFunctions.Carlson_rj(x, y, z, p), ref, tol, "Carlson_rj(" .. x .. "," .. y .. "," .. z .. "," .. p .. ")")
end

-- M5: RJ invariant under permutation of its first three arguments.
do
    local x, y, z, p = 2, 5, 11, 3
    local base = SpecialFunctions.Carlson_rj(x, y, z, p)
    local perms = {{x, z, y}, {y, x, z}, {y, z, x}, {z, x, y}, {z, y, x}}
    for _, perm in ipairs(perms) do
        Testing.Assert_close(SpecialFunctions.Carlson_rj(perm[1], perm[2], perm[3], p), base, 1e-9,
            "Carlson_rj is permutation-invariant in its first three arguments")
    end
end

-- Equal-argument closed form: RJ(a,a,a,a) = a^-1.5.
for _, a in ipairs({0.5, 1, 3.7, 20}) do
    Testing.Assert_close(SpecialFunctions.Carlson_rj(a, a, a, a), a ^ -1.5, 1e-9, "Carlson_rj(" .. a .. ",...) = a^-1.5")
end

-- M8: RD(x,y,z) = RJ(x,y,z,z) -- the identity the planned
-- E(phi|m)/E(m) implementations will lean on directly (no separate RD
-- kernel is implemented; RJ(x,y,z,z) IS RD(x,y,z)). Checked here
-- against mpmath's own elliprd as an independent reference.
local carlson_rd_via_rj_ref = {
    {2, 5, 11, 0.052525589664171590}, {0.5, 3, 8, 0.10685539016330534},
}
for _, case in ipairs(carlson_rd_via_rj_ref) do
    local x, y, z, ref = case[1], case[2], case[3], case[4]
    Testing.Assert_close(SpecialFunctions.Carlson_rj(x, y, z, z), ref, math.abs(ref) * 1e-9,
        "Carlson_rj(x,y,z,z) == RD(x,y,z) at x=" .. x .. ",y=" .. y .. ",z=" .. z)
end

-- Regression: extreme-magnitude arguments (equal in all four slots, so
-- the true value is representable only via the a^-1.5 identity itself)
-- must not fail with a spurious division-by-zero deep inside the
-- duplication loop -- caught directly in a Python model of the
-- unscaled core before this file's power-of-two scaling wrapper was
-- added.
Testing.Assert_close(SpecialFunctions.Carlson_rj(1e100, 1e100, 1e100, 1e100), (1e100) ^ -1.5, (1e100) ^ -1.5 * 1e-9,
    "Carlson_rj stays accurate (and does not fail) at extreme equal magnitude")
Testing.Assert_close(SpecialFunctions.Carlson_rf(1e150, 1e150, 1e150), (1e150) ^ -0.5, (1e150) ^ -0.5 * 1e-9,
    "Carlson_rf stays accurate at extreme magnitude")

Testing.Assert_error(function() return SpecialFunctions.Carlson_rj(-1, 1, 1, 1) end, "Carlson_rj: x, y, and z must all be nonnegative")
Testing.Assert_error(function() return SpecialFunctions.Carlson_rj(0, 0, 1, 1) end, "Carlson_rj: at most one of x, y, z may be zero")
Testing.Assert_error(function() return SpecialFunctions.Carlson_rj(1, 1, 1, 0) end, "Carlson_rj: p must be positive")
Testing.Assert_error(function() return SpecialFunctions.Carlson_rj(1, 1, 1, -1) end, "Carlson_rj: p must be positive")
Testing.Assert_error(function() return SpecialFunctions.Carlson_rj("1", 1, 1, 1) end, "Carlson_rj: x")

-- ===== Smoke test: the full third-kind reduction from the reference
-- document (test vector T20, stress-testing large negative n and m
-- together, "Cancellation region" per the doc's own randomized-test
-- strata) -- exercises all three kernels together ahead of Pi(n;phi|m)
-- itself being implemented. =====
do
    local n, phi, m = -1e10, 1e-4, -1e8
    local s, c = math.sin(phi), math.cos(phi)
    local x, y, z, p = c * c, 1 - m * s * s, 1, 1 - n * s * s
    local rf = SpecialFunctions.Carlson_rf(x, y, z)
    local rj = SpecialFunctions.Carlson_rj(x, y, z, p)
    local pi_value = s * rf + (n * s ^ 3 / 3) * rj
    Testing.Assert_close(pi_value, 1.436810311480912274254415e-5, 1e-9, "T20 smoke test (Pi(n;phi|m) via RF+RJ)")
end

-- ===== Elliptic_f / Elliptic_e_incomplete =====
--
-- Reference values below are mpmath's ellipf()/ellipe() at 30-digit
-- working precision (an independent, trusted reference), not values
-- this implementation itself produced.

local elliptic_f_ref = {
    {0.7, 0.3, 0.71651771598539313}, {1.2, 0.8, 1.4884956889493300}, {0.5, 0.5, 0.51046713562800476},
    {1.0, -2.0, 0.82956088578834132}, {1.4, 0.99, 2.3875754616023318}, {-0.9, 0.4, -0.94699792240401120},
    {0, 0.5, 0.0}, {math.pi / 2, 0.6, 1.9495677498060259},
}
for _, case in ipairs(elliptic_f_ref) do
    local phi, m, ref = case[1], case[2], case[3]
    local tol = math.max(math.abs(ref) * 1e-9, 1e-9)
    Testing.Assert_close(SpecialFunctions.Elliptic_f(phi, m), ref, tol, "Elliptic_f(" .. phi .. "," .. m .. ")")
end

local elliptic_e_incomplete_ref = {
    {0.7, 0.3, 0.68414060780670029}, {1.2, 0.8, 0.99887463983842525}, {0.5, 0.5, 0.48991095979251716},
    {1.0, -2.0, 1.2303948166988862}, {1.4, 0.99, 0.99264733285385111}, {-0.9, 0.4, -0.85692235700701177},
    {0, 0.5, 0.0}, {math.pi / 2, 0.6, 1.2984280350469132},
}
for _, case in ipairs(elliptic_e_incomplete_ref) do
    local phi, m, ref = case[1], case[2], case[3]
    local tol = math.max(math.abs(ref) * 1e-9, 1e-9)
    Testing.Assert_close(SpecialFunctions.Elliptic_e_incomplete(phi, m), ref, tol, "Elliptic_e_incomplete(" .. phi .. "," .. m .. ")")
end

-- Elliptic_f(-phi,m) = -Elliptic_f(phi,m); same for Elliptic_e_incomplete.
for _, case in ipairs({{0.8, 0.3}, {1.1, -1.5}}) do
    local phi, m = case[1], case[2]
    Testing.Assert_close(SpecialFunctions.Elliptic_f(-phi, m), -SpecialFunctions.Elliptic_f(phi, m), 1e-12,
        "Elliptic_f(-phi,m) = -Elliptic_f(phi,m)")
    Testing.Assert_close(SpecialFunctions.Elliptic_e_incomplete(-phi, m), -SpecialFunctions.Elliptic_e_incomplete(phi, m), 1e-12,
        "Elliptic_e_incomplete(-phi,m) = -Elliptic_e_incomplete(phi,m)")
end

Testing.Assert_error(function() return SpecialFunctions.Elliptic_f(2, 0.5) end, "Elliptic_f: phi must be in")
Testing.Assert_error(function() return SpecialFunctions.Elliptic_f(0.5, 1) end, "Elliptic_f: m must be < 1")
Testing.Assert_error(function() return SpecialFunctions.Elliptic_f(0.5, 1.5) end, "Elliptic_f: m must be < 1")
Testing.Assert_error(function() return SpecialFunctions.Elliptic_f("1", 0.5) end, "Elliptic_f: phi")
Testing.Assert_error(function() return SpecialFunctions.Elliptic_e_incomplete(2, 0.5) end, "Elliptic_e_incomplete: phi must be in")
Testing.Assert_error(function() return SpecialFunctions.Elliptic_e_incomplete(0.5, 1) end, "Elliptic_e_incomplete: m must be < 1")

-- ===== Elliptic_k / Elliptic_e =====
--
-- Reference values below are mpmath's ellipk()/ellipe() at 30-digit
-- working precision.

local elliptic_k_ref = {
    {0.3, 1.7138894481787911}, {0.8, 2.2572053268208538}, {-2.0, 1.1714200841467699},
    {0.99, 3.6956373629898742}, {-100, 0.36821924860914103}, {0, 1.5707963267948966},
}
for _, case in ipairs(elliptic_k_ref) do
    local m, ref = case[1], case[2]
    Testing.Assert_close(SpecialFunctions.Elliptic_k(m), ref, math.abs(ref) * 1e-9, "Elliptic_k(" .. m .. ")")
end

local elliptic_e_ref = {
    {0.3, 1.4453630644126653}, {0.8, 1.1784899243278385}, {-2.0, 2.1844381427462012},
    {0.99, 1.0159935450252239}, {-100, 10.209260919814572}, {0, 1.5707963267948966},
}
for _, case in ipairs(elliptic_e_ref) do
    local m, ref = case[1], case[2]
    Testing.Assert_close(SpecialFunctions.Elliptic_e(m), ref, math.abs(ref) * 1e-9, "Elliptic_e(" .. m .. ")")
end

-- Elliptic_k(0) = Elliptic_e(0) = pi/2 (the m=0 elementary closed form).
Testing.Assert_close(SpecialFunctions.Elliptic_k(0), math.pi / 2, 1e-12, "Elliptic_k(0) = pi/2")
Testing.Assert_close(SpecialFunctions.Elliptic_e(0), math.pi / 2, 1e-12, "Elliptic_e(0) = pi/2")

-- Consistency: Elliptic_f/Elliptic_e_incomplete at phi=pi/2 must equal
-- the complete Elliptic_k/Elliptic_e (the reference document's own
-- definition, eq. 1.2 -- the complete integral IS the incomplete one
-- evaluated at the endpoint).
for _, m in ipairs({0.3, -2, 0.99, -100}) do
    Testing.Assert_close(SpecialFunctions.Elliptic_f(math.pi / 2, m), SpecialFunctions.Elliptic_k(m), 1e-9,
        "Elliptic_f(pi/2,m) == Elliptic_k(m) at m=" .. m)
    Testing.Assert_close(SpecialFunctions.Elliptic_e_incomplete(math.pi / 2, m), SpecialFunctions.Elliptic_e(m), 1e-9,
        "Elliptic_e_incomplete(pi/2,m) == Elliptic_e(m) at m=" .. m)
end

Testing.Assert_error(function() return SpecialFunctions.Elliptic_k(1) end, "Elliptic_k: m must be < 1")
Testing.Assert_error(function() return SpecialFunctions.Elliptic_k(1.5) end, "Elliptic_k: m must be < 1")
Testing.Assert_error(function() return SpecialFunctions.Elliptic_k("1") end, "Elliptic_k: m")
Testing.Assert_error(function() return SpecialFunctions.Elliptic_e(1) end, "Elliptic_e: m must be < 1")
