-- [Discrete Math]
-- Numeral-base conversion, GCD/LCM, and primality testing -- the general
-- discrete-math primitives that don't belong to any one narrower domain
-- (combinatorics-specific functions like Pascal's triangle, binomial
-- coefficients, and Fibonacci live in the sibling GABAS-COMBINATORICS
-- project instead -- see that module's own header for why they were
-- split out rather than kept here).
--
-- Claude: not yet split into its own set of sub-submodules (base
-- conversion vs. number theory) -- deliberately, per Nacho's own stated
-- plan: get the functions working and verified first, decide the finer
-- internal split later, the same "monolith first, submodules once it
-- earns them" path GABAS_CALC_ONE_VAR.lua itself took.
local Core = require("GABAS-DISCRETE-MATH.modules.core")

local DIGIT_CHARS = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"
local DIGIT_VALUE = {}
for i = 1, #DIGIT_CHARS do
    DIGIT_VALUE[DIGIT_CHARS:sub(i, i)] = i - 1
end

-- To_base(n, base) -> string
--
-- The base-`base` representation of the integer n (any base from 2 to
-- 36, digits 0-9 then A-Z) -- e.g. To_base(202, 12) is the duodecimal
-- "14A". Via the classical repeated-division algorithm: each successive
-- remainder is the next digit, LEAST significant first, so the digits
-- collected are reversed once division reaches 0. A negative n is
-- handled as a sign followed by the representation of |n|.
--
-- Claude: this ONE verified algorithm is what every named conversion
-- wrapper below (Decimal_to_binary, Decimal_to_hexadecimal, ...) reduces
-- to -- deliberately, rather than six independently hand-written copies
-- of the same repeated-division logic (see the header note on
-- Binary_to_hexadecimal below for the same reasoning applied to
-- base-to-base conversions).
local function To_base(n, base)
    Core.assert_integer(n, "To_base: n")
    Core.assert_valid_base(base, "To_base: base")
    if n == 0 then
        return "0"
    end
    local negative = n < 0
    local magnitude = negative and -n or n
    local digits = {}
    while magnitude > 0 do
        local remainder = magnitude % base
        digits[#digits + 1] = DIGIT_CHARS:sub(remainder + 1, remainder + 1)
        magnitude = magnitude // base
    end
    local out = {}
    for i = #digits, 1, -1 do
        out[#out + 1] = digits[i]
    end
    local result = table.concat(out)
    return negative and ("-" .. result) or result
end

-- From_base(str, base) -> integer
--
-- The inverse of To_base: parses a string of digits (0-9, then A-Z,
-- case-insensitive) in the given base back into a plain integer, via
-- Horner's rule (each new digit shifts the value accumulated so far up
-- by one base-place and adds the digit). An optional leading "-" is
-- accepted for a negative value.
--
-- Claude: every character is validated to be an actual digit BELOW
-- `base` (e.g. "2" is not a valid binary digit) BEFORE any arithmetic
-- happens, so a malformed string is rejected with the exact offending
-- character named, not silently misparsed or left to a confusing
-- downstream failure. The accumulation is also checked for 64-bit
-- integer overflow at every step, since Lua 5.4's integer type wraps
-- around silently rather than erroring on overflow -- exactly the
-- silent-wrong-answer failure mode this whole project refuses to allow.
local function From_base(str, base)
    Core.assert_nonblank_string(str, "From_base: str")
    Core.assert_valid_base(base, "From_base: base")
    local negative = false
    local start = 1
    if str:sub(1, 1) == "-" then
        negative = true
        start = 2
    end
    assert(start <= #str, "From_base: \"" .. str .. "\" has no digits after the sign.")
    local value = 0
    for i = start, #str do
        local c = str:sub(i, i):upper()
        local digit = DIGIT_VALUE[c]
        assert(digit ~= nil, "From_base: \"" .. str .. "\" contains \"" .. str:sub(i, i) ..
            "\", which is not a valid digit character.")
        assert(digit < base, "From_base: \"" .. str .. "\" contains the digit \"" .. str:sub(i, i) ..
            "\", which is not valid in base " .. base .. ".")
        assert(value <= (math.maxinteger - digit) // base,
            "From_base: \"" .. str .. "\" overflows a 64-bit integer.")
        value = value * base + digit
    end
    return negative and -value or value
end

-- Claude: each of these is a thin, explicitly-named wrapper around
-- To_base/From_base -- one verified algorithm underneath, not several
-- independently-maintained copies of the same logic (the exact
-- duplication risk Nacho flagged when we discarded putting a function in
-- two modules at once).
local function Decimal_to_binary(n) return To_base(n, 2) end
local function Binary_to_decimal(str) return From_base(str, 2) end
local function Decimal_to_octal(n) return To_base(n, 8) end
local function Octal_to_decimal(str) return From_base(str, 8) end
local function Decimal_to_hexadecimal(n) return To_base(n, 16) end
local function Hexadecimal_to_decimal(str) return From_base(str, 16) end
local function Decimal_to_duodecimal(n) return To_base(n, 12) end
local function Duodecimal_to_decimal(str) return From_base(str, 12) end

-- Claude: a direct base-to-base conversion goes THROUGH the same
-- verified To_base/From_base pair (parse in the source base, re-render
-- in the target base) rather than a separate bit-grouping shortcut
-- (e.g. exploiting that 16 = 2^4). At the integer sizes this module
-- deals with, there is no meaningful performance reason to maintain a
-- second, different algorithm that would need its own separate
-- verification just to save one intermediate value.
local function Binary_to_hexadecimal(str) return To_base(From_base(str, 2), 16) end
local function Hexadecimal_to_binary(str) return To_base(From_base(str, 16), 2) end

-- Gcd(a, b) -> the greatest common divisor of a and b, via the Euclidean
-- algorithm (repeated remainder, ~2500 years old and about as settled as
-- an algorithm gets). Always returns a nonnegative value; Gcd(0, 0) is 0
-- by convention (every integer divides 0, so there is no single greatest
-- one unless this base case is fixed).
local function Gcd(a, b)
    Core.assert_integer(a, "Gcd: a")
    Core.assert_integer(b, "Gcd: b")
    a, b = math.abs(a), math.abs(b)
    while b ~= 0 do
        a, b = b, a % b
    end
    return a
end

-- Lcm(a, b) -> the least common multiple of a and b, via
-- (a/gcd(a,b))*b -- dividing by the gcd FIRST, rather than computing
-- a*b and dividing after, keeps the intermediate product as small as
-- possible, the standard defense against unnecessary overflow. Lcm(0,
-- anything) is 0 by convention.
local function Lcm(a, b)
    Core.assert_integer(a, "Lcm: a")
    Core.assert_integer(b, "Lcm: b")
    if a == 0 or b == 0 then
        return 0
    end
    local scaled = math.abs(a) // Gcd(a, b)
    assert(scaled <= math.maxinteger // math.abs(b),
        "Lcm: the result overflows a 64-bit integer.")
    return scaled * math.abs(b)
end

-- Is_prime(n) -> true iff n is a prime number, via trial division up to
-- sqrt(n) (skipping even candidates after 2), the standard elementary
-- primality test -- adequate for this module's scope. An arbitrary-
-- precision or probabilistic test (e.g. Miller-Rabin, for
-- cryptographic-scale numbers) is a different, much larger undertaking,
-- not needed here.
local function Is_prime(n)
    Core.assert_integer(n, "Is_prime: n")
    if n < 2 then return false end
    if n == 2 or n == 3 then return true end
    if n % 2 == 0 then return false end
    local i = 3
    while i * i <= n do
        if n % i == 0 then return false end
        i = i + 2
    end
    return true
end

return {
    To_base = To_base,
    From_base = From_base,
    Decimal_to_binary = Decimal_to_binary,
    Binary_to_decimal = Binary_to_decimal,
    Decimal_to_octal = Decimal_to_octal,
    Octal_to_decimal = Octal_to_decimal,
    Decimal_to_hexadecimal = Decimal_to_hexadecimal,
    Hexadecimal_to_decimal = Hexadecimal_to_decimal,
    Decimal_to_duodecimal = Decimal_to_duodecimal,
    Duodecimal_to_decimal = Duodecimal_to_decimal,
    Binary_to_hexadecimal = Binary_to_hexadecimal,
    Hexadecimal_to_binary = Hexadecimal_to_binary,
    Gcd = Gcd,
    Lcm = Lcm,
    Is_prime = Is_prime,
}
