local SpecialFunctions = require("GABAS-SPECIAL-FUNCTIONS.GABAS_SPECIAL_FUNCTIONS")
local Testing = require("GABAS-LINAL.modules.testing")

-- ===== Chebyshev_t =====
--
-- Reference values below come from the independent closed forms
-- T_n(x) = cos(n*arccos(x)) for |x|<=1 and T_n(x) = cosh(n*arccosh(x))
-- for x>1, computed via mpmath at 30-digit precision -- not from this
-- implementation's own recurrence, and not from mpmath's own chebyt
-- (whose internal hypergeometric representation fails to converge for
-- some of these inputs -- an mpmath quirk, not a property of T_n
-- itself, worth noting since it's the kind of thing that could
-- otherwise look like a bug in this implementation during review).

local chebyshev_t_ref = {
    {0, 0.5, 1.0}, {1, 0.5, 0.5}, {2, 0.5, -0.5}, {3, 0.5, -1.0},
    {5, 0.3, 0.99888}, {10, 0.9, -0.2007474688000005},
    {0, 2, 1.0}, {1, 2, 2.0}, {3, 2, 26.0}, {5, 1.5, 61.5},
    {10, 5, 4517251249.0}, {20, 1.2, 127274.68371953907},
}
for _, case in ipairs(chebyshev_t_ref) do
    local n, x, ref = case[1], case[2], case[3]
    local tol = math.max(math.abs(ref) * 1e-9, 1e-9)
    Testing.Assert_close(SpecialFunctions.Chebyshev_t(n, x), ref, tol, "Chebyshev_t(" .. n .. "," .. x .. ")")
end

-- Closed forms at the domain endpoints: T_n(1) = 1, T_n(-1) = (-1)^n.
for n = 0, 10 do
    assert(SpecialFunctions.Chebyshev_t(n, 1) == 1, "Chebyshev_t(" .. n .. ",1) == 1")
    local expected = (n % 2 == 0) and 1 or -1
    assert(SpecialFunctions.Chebyshev_t(n, -1) == expected, "Chebyshev_t(" .. n .. ",-1) == (-1)^" .. n)
end

-- T_n(-x) = (-1)^n * T_n(x).
for _, n in ipairs({0, 1, 2, 3, 7}) do
    local x = 1.7
    local expected = ((n % 2 == 0) and 1 or -1) * SpecialFunctions.Chebyshev_t(n, x)
    Testing.Assert_close(SpecialFunctions.Chebyshev_t(n, -x), expected, 1e-9,
        "Chebyshev_t(" .. n .. ",-x) = (-1)^" .. n .. "*Chebyshev_t(" .. n .. ",x)")
end

-- Regression: an integer-looking argument (no decimal point) is a Lua
-- *integer*, and without forcing float arithmetic up front, this
-- recurrence would silently wrap around 64-bit integer overflow
-- instead of promoting to float once |T_n(x)| exceeds ~9.2e18 --
-- Chebyshev_t(34, -2) (called with plain integer literals, exactly as
-- a user would naturally write it) came back negative and wildly wrong
-- under that version.
Testing.Assert_close(SpecialFunctions.Chebyshev_t(34, -2), 13969685227624439047.0, 13969685227624439047.0 * 1e-9,
    "Chebyshev_t(34,-2) does not wrap from integer overflow")

Testing.Assert_error(function() return SpecialFunctions.Chebyshev_t(-1, 0.5) end, "Chebyshev_t: n must be a nonnegative integer")
Testing.Assert_error(function() return SpecialFunctions.Chebyshev_t(1.5, 0.5) end, "Chebyshev_t: n must be a nonnegative integer")
Testing.Assert_error(function() return SpecialFunctions.Chebyshev_t(0, "1") end, "Chebyshev_t: x")

-- ===== Chebyshev_u =====
--
-- Reference values below come from the independent closed forms
-- U_n(x) = sin((n+1)*arccos(x))/sin(arccos(x)) for |x|<1 and the
-- analogous sinh/arccosh form for |x|>1, computed via mpmath at
-- 30-digit precision, for the same reason as Chebyshev_t's above.

local chebyshev_u_ref = {
    {0, 0.5, 1.0}, {1, 0.5, 1.0}, {2, 0.5, 0.0}, {3, 0.5, -1.0},
    {5, 0.3, 1.01376}, {10, 0.9, -2.2234571776000006},
    {0, 2, 1.0}, {1, 2, 4.0}, {3, 2, 56.0}, {5, 1.5, 144.0},
    {10, 5, 9127651499.0}, {20, 1.2, 357523.24982635909},
}
for _, case in ipairs(chebyshev_u_ref) do
    local n, x, ref = case[1], case[2], case[3]
    local tol = math.max(math.abs(ref) * 1e-9, 1e-9)
    Testing.Assert_close(SpecialFunctions.Chebyshev_u(n, x), ref, tol, "Chebyshev_u(" .. n .. "," .. x .. ")")
end

-- Closed forms at the domain endpoints: U_n(1) = n+1, U_n(-1) = (-1)^n*(n+1).
for n = 0, 10 do
    assert(SpecialFunctions.Chebyshev_u(n, 1) == n + 1, "Chebyshev_u(" .. n .. ",1) == " .. (n + 1))
    local expected = ((n % 2 == 0) and 1 or -1) * (n + 1)
    assert(SpecialFunctions.Chebyshev_u(n, -1) == expected, "Chebyshev_u(" .. n .. ",-1) == (-1)^" .. n .. "*(n+1)")
end

-- U_n(-x) = (-1)^n * U_n(x).
for _, n in ipairs({0, 1, 2, 3, 7}) do
    local x = 1.7
    local expected = ((n % 2 == 0) and 1 or -1) * SpecialFunctions.Chebyshev_u(n, x)
    Testing.Assert_close(SpecialFunctions.Chebyshev_u(n, -x), expected, 1e-6,
        "Chebyshev_u(" .. n .. ",-x) = (-1)^" .. n .. "*Chebyshev_u(" .. n .. ",x)")
end

-- Regression: same integer-overflow hazard as Chebyshev_t above.
Testing.Assert_close(SpecialFunctions.Chebyshev_u(34, -2), 30100488280951055759.0, 30100488280951055759.0 * 1e-9,
    "Chebyshev_u(34,-2) does not wrap from integer overflow")

Testing.Assert_error(function() return SpecialFunctions.Chebyshev_u(-1, 0.5) end, "Chebyshev_u: n must be a nonnegative integer")
Testing.Assert_error(function() return SpecialFunctions.Chebyshev_u(1.5, 0.5) end, "Chebyshev_u: n must be a nonnegative integer")
Testing.Assert_error(function() return SpecialFunctions.Chebyshev_u(0, "1") end, "Chebyshev_u: x")

-- ===== Legendre_p =====
--
-- Reference values below are mpmath's legendre() at 30-digit working
-- precision (an independent, trusted reference), not values this
-- implementation itself produced.

local legendre_p_ref = {
    {0, 0.5, 1.0}, {1, 0.5, 0.5}, {2, 0.5, -0.125}, {3, 0.5, -0.4375},
    {5, 0.3, 0.34538625}, {10, 0.9, -0.26314561785585953},
    {0, 2, 1.0}, {1, 2, 2.0}, {3, 2, 17.0}, {5, 1.5, 33.08203125},
    {10, 5, 1600472677.0}, {20, 1.2, 38022.330892378732},
}
for _, case in ipairs(legendre_p_ref) do
    local n, x, ref = case[1], case[2], case[3]
    local tol = math.max(math.abs(ref) * 1e-9, 1e-9)
    Testing.Assert_close(SpecialFunctions.Legendre_p(n, x), ref, tol, "Legendre_p(" .. n .. "," .. x .. ")")
end

-- Closed forms at the domain endpoints: P_n(1) = 1, P_n(-1) = (-1)^n.
for n = 0, 10 do
    assert(SpecialFunctions.Legendre_p(n, 1) == 1, "Legendre_p(" .. n .. ",1) == 1")
    local expected = (n % 2 == 0) and 1 or -1
    assert(SpecialFunctions.Legendre_p(n, -1) == expected, "Legendre_p(" .. n .. ",-1) == (-1)^" .. n)
end

-- P_n(-x) = (-1)^n * P_n(x).
for _, n in ipairs({0, 1, 2, 3, 7}) do
    local x = 1.7
    local expected = ((n % 2 == 0) and 1 or -1) * SpecialFunctions.Legendre_p(n, x)
    Testing.Assert_close(SpecialFunctions.Legendre_p(n, -x), expected, 1e-9,
        "Legendre_p(" .. n .. ",-x) = (-1)^" .. n .. "*Legendre_p(" .. n .. ",x)")
end

-- Same integer-only-argument case that exposed the Chebyshev_t/
-- Chebyshev_u overflow bug -- Legendre_p's own division-based recurrence
-- self-promotes to float regardless (see its header), so this is a
-- spot check confirming that holds, not a regression for a bug that
-- was actually caught here.
Testing.Assert_close(SpecialFunctions.Legendre_p(34, -2), 2797276292689877223.3, 2797276292689877223.3 * 1e-9,
    "Legendre_p(34,-2) stays accurate with integer-only arguments")

Testing.Assert_error(function() return SpecialFunctions.Legendre_p(-1, 0.5) end, "Legendre_p: n must be a nonnegative integer")
Testing.Assert_error(function() return SpecialFunctions.Legendre_p(1.5, 0.5) end, "Legendre_p: n must be a nonnegative integer")
Testing.Assert_error(function() return SpecialFunctions.Legendre_p(0, "1") end, "Legendre_p: x")

-- ===== Hermite_h =====
--
-- Reference values below are mpmath's hermite() (physicists'
-- convention) at 30-digit working precision (an independent, trusted
-- reference), not values this implementation itself produced.

local hermite_h_ref = {
    {0, 0.5, 1.0}, {1, 0.5, 1.0}, {2, 0.5, -1.0}, {3, 0.5, -5.0},
    {5, 0.3, 31.75776}, {10, 0.9, 26314.366684262397},
    {0, 2, 1.0}, {1, 2, 4.0}, {3, 2, 40.0}, {5, 1.5, -117.0},
    {6, 5, 717880.0}, {8, 1.2, 594.57048575999936},
}
for _, case in ipairs(hermite_h_ref) do
    local n, x, ref = case[1], case[2], case[3]
    local tol = math.max(math.abs(ref) * 1e-9, 1e-9)
    Testing.Assert_close(SpecialFunctions.Hermite_h(n, x), ref, tol, "Hermite_h(" .. n .. "," .. x .. ")")
end

-- Closed forms at x=0: H_n(0)=0 for odd n; a few known even-n values.
for _, n in ipairs({1, 3, 5, 7}) do
    assert(SpecialFunctions.Hermite_h(n, 0) == 0, "Hermite_h(" .. n .. ",0) == 0")
end
assert(SpecialFunctions.Hermite_h(0, 0) == 1, "Hermite_h(0,0) == 1")
assert(SpecialFunctions.Hermite_h(2, 0) == -2, "Hermite_h(2,0) == -2")
assert(SpecialFunctions.Hermite_h(4, 0) == 12, "Hermite_h(4,0) == 12")

-- H_n(-x) = (-1)^n * H_n(x).
for _, n in ipairs({0, 1, 2, 3, 7}) do
    local x = 1.7
    local expected = ((n % 2 == 0) and 1 or -1) * SpecialFunctions.Hermite_h(n, x)
    Testing.Assert_close(SpecialFunctions.Hermite_h(n, -x), expected, 1e-6,
        "Hermite_h(" .. n .. ",-x) = (-1)^" .. n .. "*Hermite_h(" .. n .. ",x)")
end

-- Same integer-only-argument spot check as Legendre_p above.
Testing.Assert_close(SpecialFunctions.Hermite_h(34, -2), 4.5791990566263883069e+24, 4.5791990566263883069e+24 * 1e-9,
    "Hermite_h(34,-2) stays accurate with integer-only arguments")

Testing.Assert_error(function() return SpecialFunctions.Hermite_h(-1, 0.5) end, "Hermite_h: n must be a nonnegative integer")
Testing.Assert_error(function() return SpecialFunctions.Hermite_h(1.5, 0.5) end, "Hermite_h: n must be a nonnegative integer")
Testing.Assert_error(function() return SpecialFunctions.Hermite_h(0, "1") end, "Hermite_h: x")
