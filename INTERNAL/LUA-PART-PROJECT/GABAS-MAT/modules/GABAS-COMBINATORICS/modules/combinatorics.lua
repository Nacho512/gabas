-- [Combinatorics]
-- Factorial, binomial coefficients, Pascal's triangle, and Fibonacci --
-- the specifically combinatorial/recurrence-based discrete-math
-- primitives, split out from the more general GABAS-DISCRETE-MATH
-- project (numeral-base conversion, GCD/LCM, primality) because they
-- belong to a narrower, more specific domain (counting and recurrences,
-- not discrete math broadly) -- see that project's own header for the
-- other half of this split.
--
-- Claude: not yet split into its own set of sub-submodules -- same
-- "monolith first, submodules once it earns them" plan as
-- GABAS-DISCRETE-MATH.
--
-- Every function below guards against silent 64-bit integer overflow --
-- Lua 5.4's integer type wraps around rather than erroring when a
-- computation exceeds its range, and factorials/binomial
-- coefficients/Fibonacci numbers all grow fast enough to hit that
-- ceiling for very ordinary-looking inputs (20! is already the largest
-- representable factorial). A caller must get an honest error there,
-- never a silently wrong wrapped-around number.
local Core = require("GABAS-COMBINATORICS.modules.core")

-- Factorial(n) -> n! (0! = 1 by convention), for a nonnegative integer n,
-- via a plain iterative product -- no recursion, so there is no stack-
-- depth concern for larger n. The largest representable factorial in a
-- signed 64-bit integer is 20! (2432902008176640000); 21! would already
-- overflow, so this is checked before every multiplication, not
-- discovered after the fact.
local function Factorial(n)
    Core.assert_nonneg_integer(n, "Factorial: n")
    local result = 1
    for i = 2, n do
        assert(result <= math.maxinteger // i,
            "Factorial: " .. n .. "! overflows a 64-bit integer (the largest representable factorial is 20!).")
        result = result * i
    end
    return result
end

-- Binomial_coefficient(n, k) -> C(n,k), the number of ways to choose k
-- items from n, for a nonnegative integer n and any integer k.
-- Conventionally 0 whenever k < 0 or k > n -- not an error; this is the
-- standard convention that lets binomial sums and identities work
-- without special-casing their bounds.
--
-- Claude: computed via the classical multiplicative formula --
-- C(n,k) = product_(i=1..k) (n-k+i)/i -- rather than n!/(k!(n-k)!),
-- which would overflow far sooner (this module's own Factorial already
-- tops out at 20!, but C(n,k) itself stays representable for much
-- larger n, as long as the FINAL result fits). Always iterates over
-- k' = min(k, n-k) (the C(n,k) = C(n,n-k) symmetry), halving the work
-- when k is close to n. Each partial product is guaranteed to be an
-- exact integer at every step -- a classical number-theory fact (the
-- product of any k consecutive integers is divisible by k!) -- so the
-- integer division by i below never discards a fractional part; this is
-- verified, not just assumed, by the exhaustive small-(n,k) cross-check
-- against Pascal's triangle in the test suite.
local function Binomial_coefficient(n, k)
    Core.assert_nonneg_integer(n, "Binomial_coefficient: n")
    Core.assert_integer(k, "Binomial_coefficient: k")
    if k < 0 or k > n then
        return 0
    end
    local kk = math.min(k, n - k)
    local result = 1
    for i = 1, kk do
        local numerator_factor = n - kk + i
        assert(result <= math.maxinteger // numerator_factor,
            "Binomial_coefficient: C(" .. n .. "," .. k .. ") overflows a 64-bit integer.")
        result = (result * numerator_factor) // i
    end
    return result
end

-- Pascal_row(n) -> a table where row[k] == Binomial_coefficient(n, k)
-- for k = 0, 1, ..., n -- Pascal's triangle's nth row.
--
-- Claude: deliberately 0-indexed (row[0], row[1], ..., row[n]) -- an
-- explicit exception to Lua's usual 1-based convention, chosen because
-- row[k] meaning literally C(n,k), with no off-by-one translation, is
-- worth the departure. Plain numeric access (row[k]) works fine; just
-- don't reach for ipairs() on it (which starts at index 1 and would
-- silently skip row[0]).
local function Pascal_row(n)
    Core.assert_nonneg_integer(n, "Pascal_row: n")
    local row = {}
    for k = 0, n do
        row[k] = Binomial_coefficient(n, k)
    end
    return row
end

-- Pascal_triangle(max_row) -> a table where triangle[i] == Pascal_row(i)
-- for i = 0, 1, ..., max_row -- 0-indexed the same way, and for the same
-- reason, as Pascal_row itself.
local function Pascal_triangle(max_row)
    Core.assert_nonneg_integer(max_row, "Pascal_triangle: max_row")
    local triangle = {}
    for i = 0, max_row do
        triangle[i] = Pascal_row(i)
    end
    return triangle
end

-- Fibonacci(n) -> the nth Fibonacci number (F(0) = 0, F(1) = 1,
-- F(n) = F(n-1) + F(n-2)), computed iteratively (O(n) additions, no
-- recursion) -- more than sufficient for this module's current scope. A
-- fast-doubling O(log n) algorithm is a real future increment if very
-- large n is ever actually needed, not folded in prematurely here. F(92)
-- is the largest Fibonacci number that still fits in a signed 64-bit
-- integer; F(93) does not, checked before it happens rather than
-- discovered after.
local function Fibonacci(n)
    Core.assert_nonneg_integer(n, "Fibonacci: n")
    if n == 0 then return 0 end
    local a, b = 0, 1
    for _ = 2, n do
        assert(b <= math.maxinteger - a,
            "Fibonacci: F(" .. n .. ") overflows a 64-bit integer (the largest representable is F(92)).")
        a, b = b, a + b
    end
    return b
end

return {
    Factorial = Factorial,
    Binomial_coefficient = Binomial_coefficient,
    Pascal_row = Pascal_row,
    Pascal_triangle = Pascal_triangle,
    Fibonacci = Fibonacci,
}
