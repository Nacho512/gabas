local Combinatorics = require("GABAS-COMBINATORICS.GABAS_COMBINATORICS")
local Testing = require("GABAS-LINAL.modules.testing")

-- ===== Factorial =====

assert(Combinatorics.Factorial(0) == 1, "0! == 1 by convention")
assert(Combinatorics.Factorial(1) == 1, "1! == 1")
assert(Combinatorics.Factorial(5) == 120, "5! == 120")
assert(Combinatorics.Factorial(20) == 2432902008176640000, "20! is the largest representable factorial")

Testing.Assert_error(function() return Combinatorics.Factorial(21) end, "overflows a 64-bit integer")
Testing.Assert_error(function() return Combinatorics.Factorial(-1) end, "Factorial: n")
Testing.Assert_error(function() return Combinatorics.Factorial(2.5) end, "Factorial: n")

-- ===== Binomial_coefficient =====
--
-- Reference values below (C(52,5), C(20,10)) are independent, well-known
-- closed-form results (5-card poker hands out of a standard deck; a
-- central binomial coefficient), not values this implementation itself
-- produced.

assert(Combinatorics.Binomial_coefficient(52, 5) == 2598960, "C(52,5), 5-card poker hands")
assert(Combinatorics.Binomial_coefficient(20, 10) == 184756, "C(20,10)")
assert(Combinatorics.Binomial_coefficient(10, 3) == 120, "C(10,3)")
assert(Combinatorics.Binomial_coefficient(5, 0) == 1, "C(n,0) == 1")
assert(Combinatorics.Binomial_coefficient(5, 5) == 1, "C(n,n) == 1")
assert(Combinatorics.Binomial_coefficient(0, 0) == 1, "C(0,0) == 1")

-- Out-of-range k is 0 BY CONVENTION, not an error.
assert(Combinatorics.Binomial_coefficient(5, 6) == 0, "C(n,k) == 0 for k > n")
assert(Combinatorics.Binomial_coefficient(5, -1) == 0, "C(n,k) == 0 for k < 0")

-- Symmetry: C(n,k) == C(n,n-k).
for n = 0, 15 do
    for k = 0, n do
        assert(Combinatorics.Binomial_coefficient(n, k) == Combinatorics.Binomial_coefficient(n, n - k),
            "C(" .. n .. "," .. k .. ") == C(n,n-k) symmetry")
    end
end

-- Genuinely astronomical (C(100,50) has 30 digits) -- must report an
-- honest overflow error, never a silently wrapped-around 64-bit value.
Testing.Assert_error(function() return Combinatorics.Binomial_coefficient(100, 50) end,
    "overflows a 64-bit integer")

Testing.Assert_error(function() return Combinatorics.Binomial_coefficient(-1, 0) end,
    "Binomial_coefficient: n")

-- ===== Pascal_row / Pascal_triangle =====

local row5 = Combinatorics.Pascal_row(5)
for k = 0, 5 do
    assert(row5[k] == Combinatorics.Binomial_coefficient(5, k),
        "Pascal_row(5)[" .. k .. "] matches Binomial_coefficient")
end
assert(row5[0] == 1 and row5[1] == 5 and row5[2] == 10 and row5[3] == 10
    and row5[4] == 5 and row5[5] == 1, "Pascal_row(5) == {1,5,10,10,5,1}")

-- Pascal's own additive recurrence: row[k] == prev_row[k-1] + prev_row[k]
-- (the "1 1" edges are the base case) -- checked directly, an
-- independent identity, not just agreement with Binomial_coefficient.
local triangle = Combinatorics.Pascal_triangle(15)
for n = 1, 15 do
    local row, prev = triangle[n], triangle[n - 1]
    assert(row[0] == 1 and row[n] == 1, "Pascal_triangle row " .. n .. " edges are 1")
    for k = 1, n - 1 do
        assert(row[k] == prev[k - 1] + prev[k],
            "Pascal's recurrence at row " .. n .. ", k=" .. k)
    end
end

Testing.Assert_error(function() return Combinatorics.Pascal_row(-1) end, "Pascal_row: n")

-- ===== Fibonacci =====
--
-- F(0..10) below is the textbook sequence, independent of this
-- implementation.

local expected = {0, 1, 1, 2, 3, 5, 8, 13, 21, 34, 55}
for n = 0, 10 do
    assert(Combinatorics.Fibonacci(n) == expected[n + 1], "Fibonacci(" .. n .. ")")
end
assert(Combinatorics.Fibonacci(92) == 7540113804746346429, "F(92) is the largest representable Fibonacci number")

Testing.Assert_error(function() return Combinatorics.Fibonacci(93) end, "overflows a 64-bit integer")
Testing.Assert_error(function() return Combinatorics.Fibonacci(-1) end, "Fibonacci: n")

print("PASS test_combinatorics (internal)")
