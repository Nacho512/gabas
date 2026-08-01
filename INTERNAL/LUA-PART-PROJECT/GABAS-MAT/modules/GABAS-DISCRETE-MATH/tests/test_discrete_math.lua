local DiscreteMath = require("GABAS-DISCRETE-MATH.GABAS_DISCRETE_MATH")
local Testing = require("GABAS-LINAL.modules.testing")

-- ===== To_base / From_base (the general engine every named wrapper
-- below reduces to) =====

assert(DiscreteMath.To_base(202, 12) == "14A", "To_base(202,12), duodecimal")
assert(DiscreteMath.To_base(0, 16) == "0", "To_base(0, base) is always \"0\"")
assert(DiscreteMath.To_base(255, 16) == "FF", "To_base(255,16)")
assert(DiscreteMath.To_base(-255, 16) == "-FF", "To_base(-255,16), negative")
assert(DiscreteMath.To_base(35, 36) == "Z", "To_base(35,36), the largest single digit")

assert(DiscreteMath.From_base("14A", 12) == 202, "From_base(\"14A\",12)")
assert(DiscreteMath.From_base("FF", 16) == 255, "From_base(\"FF\",16)")
assert(DiscreteMath.From_base("ff", 16) == 255, "From_base is case-insensitive")
assert(DiscreteMath.From_base("-FF", 16) == -255, "From_base(\"-FF\",16), negative")
assert(DiscreteMath.From_base("0", 2) == 0, "From_base(\"0\",2)")

-- Round-trip over a spread of values and bases.
for _, base in ipairs({2, 3, 7, 8, 12, 16, 20, 36}) do
    for _, n in ipairs({0, 1, 7, 42, 1000, -1000, 123456789}) do
        assert(DiscreteMath.From_base(DiscreteMath.To_base(n, base), base) == n,
            "round-trip n=" .. n .. " base=" .. base)
    end
end

-- Invalid digit for the given base is rejected, naming the offending
-- character.
Testing.Assert_error(function() return DiscreteMath.From_base("102", 2) end,
    "not valid in base 2")
-- "G" is a genuine base-36 digit (value 16), just not a valid one IN
-- base 16 -- this exercises the "in-alphabet but out-of-range for this
-- base" branch, distinct from the next check below (not a digit at all).
Testing.Assert_error(function() return DiscreteMath.From_base("1G", 16) end,
    "not valid in base 16")
Testing.Assert_error(function() return DiscreteMath.From_base("1$", 16) end,
    "not a valid digit character")

-- Input validation.
Testing.Assert_error(function() return DiscreteMath.To_base(3.5, 10) end, "To_base: n")
Testing.Assert_error(function() return DiscreteMath.To_base(10, 1) end, "To_base: base")
Testing.Assert_error(function() return DiscreteMath.To_base(10, 37) end, "To_base: base")
Testing.Assert_error(function() return DiscreteMath.From_base("", 10) end, "From_base: str")
Testing.Assert_error(function() return DiscreteMath.From_base("-", 10) end, "no digits after the sign")

-- ===== Named conversion wrappers =====

assert(DiscreteMath.Decimal_to_binary(13) == "1101", "Decimal_to_binary(13)")
assert(DiscreteMath.Binary_to_decimal("1101") == 13, "Binary_to_decimal(\"1101\")")
assert(DiscreteMath.Decimal_to_hexadecimal(4096) == "1000", "Decimal_to_hexadecimal(4096)")
assert(DiscreteMath.Hexadecimal_to_decimal("1000") == 4096, "Hexadecimal_to_decimal(\"1000\")")
assert(DiscreteMath.Decimal_to_octal(8) == "10", "Decimal_to_octal(8)")
assert(DiscreteMath.Octal_to_decimal("10") == 8, "Octal_to_decimal(\"10\")")
assert(DiscreteMath.Decimal_to_duodecimal(144) == "100", "Decimal_to_duodecimal(144)")
assert(DiscreteMath.Duodecimal_to_decimal("100") == 144, "Duodecimal_to_decimal(\"100\")")

-- Binary <-> hexadecimal, directly.
assert(DiscreteMath.Binary_to_hexadecimal("11111111") == "FF", "Binary_to_hexadecimal(\"11111111\")")
assert(DiscreteMath.Hexadecimal_to_binary("FF") == "11111111", "Hexadecimal_to_binary(\"FF\")")

-- ===== Gcd / Lcm =====

assert(DiscreteMath.Gcd(48, 18) == 6, "Gcd(48,18)")
assert(DiscreteMath.Gcd(17, 5) == 1, "Gcd of two coprime numbers")
assert(DiscreteMath.Gcd(0, 7) == 7, "Gcd(0,n) == n")
assert(DiscreteMath.Gcd(0, 0) == 0, "Gcd(0,0) == 0 by convention")
assert(DiscreteMath.Gcd(-48, 18) == 6, "Gcd is sign-independent")

assert(DiscreteMath.Lcm(4, 6) == 12, "Lcm(4,6)")
assert(DiscreteMath.Lcm(21, 6) == 42, "Lcm(21,6)")
assert(DiscreteMath.Lcm(0, 5) == 0, "Lcm(0,n) == 0 by convention")

Testing.Assert_error(function() return DiscreteMath.Gcd(1.5, 2) end, "Gcd: a")

-- ===== Is_prime =====

for _, p in ipairs({2, 3, 5, 7, 11, 13, 97, 7919}) do
    assert(DiscreteMath.Is_prime(p), p .. " should be prime")
end
for _, c in ipairs({0, 1, 4, 6, 9, 100, 7920}) do
    assert(not DiscreteMath.Is_prime(c), c .. " should not be prime")
end
assert(not DiscreteMath.Is_prime(-7), "negative numbers are never prime")

print("PASS test_discrete_math (internal)")
