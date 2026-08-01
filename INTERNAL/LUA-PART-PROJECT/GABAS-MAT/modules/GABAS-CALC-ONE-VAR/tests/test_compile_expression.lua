local CalcOneVar = require("GABAS-CALC-ONE-VAR.GABAS_CALC_ONE_VAR")
local Testing = require("GABAS-LINAL.modules.testing")

local Compile = CalcOneVar.Compile_expression

-- ===== Identifier whitelist (Compile_Expression's own "approved syntax"
-- policy, technical document section 17.1) =====
--
-- Anything the compiled source would reference that isn't the variable x
-- or one of this package's own recognized functions/constants must be
-- rejected BEFORE load() ever runs -- a typo, a disallowed global, or a
-- Lua keyword smuggled in as if it were part of the expression.

Testing.Assert_error(function() return Compile("sni(x)") end, "Compile_Expression")
Testing.Assert_error(function() return Compile("sni(x)") end, "\"sni\"")

-- "math" was deliberately removed from the compiled environment -- a
-- user can no longer bypass the Dual/Complex-aware dispatchers by
-- reaching math.sin/math.pi/... directly.
Testing.Assert_error(function() return Compile("math.sin(x)") end, "\"math\"")
Testing.Assert_error(function() return Compile("math.pi") end, "\"math\"")
Testing.Assert_error(function() return Compile("math.random()") end, "\"math\"")

-- A completely unrecognized global function call.
Testing.Assert_error(function() return Compile("foo(x)") end, "\"foo\"")

-- Lua keywords are lexically indistinguishable from identifiers at this
-- level, and simply aren't on the whitelist -- rejected the same way.
-- (Whitespace is insignificant everywhere in this syntax -- see
-- strip_whitespace's own header -- so "x and 1" becomes "xand1" before
-- this check ever runs; the offending merged token's exact spelling
-- isn't pinned here, only that it is correctly rejected either way.)
Testing.Assert_error(function() return Compile("x and 1") end, "approved expression syntax")
Testing.Assert_error(function() return Compile("x or 0") end, "approved expression syntax")
Testing.Assert_error(function() return Compile("not(x)") end, "\"not\"")

-- ===== Numeric literals must never be misread as identifiers =====
--
-- The exponent marker in scientific notation ("e"/"E") looks exactly
-- like the start of a bare identifier -- mask_numeric_literals must
-- strip these BEFORE identifiers are scanned for, or a perfectly
-- ordinary numeral would be falsely rejected.

local f1 = Compile("1e10 * x")
Testing.Assert_close(f1(2), 2e10, 1, "1e10 * x compiles and evaluates correctly")

local f2 = Compile("x + 6.13")
Testing.Assert_close(f2(1), 7.13, 1e-9, "x + 6.13 (decimal literal) compiles correctly")

local f3 = Compile(".5e-3 + x")
Testing.Assert_close(f3(0), 0.0005, 1e-9, "leading-dot scientific notation compiles correctly")

local f4 = Compile("root_6.13(x)")
Testing.Assert_close(f4(2), 2 ^ (1 / 6.13), 1e-9, "root_6.13(x), a decimal root index, still works")

-- ===== Every recognized function/constant still works after the "math"
-- removal and the new whitelist check (no regression in the approved
-- vocabulary itself) =====

Testing.Assert_close((Compile("sin(x)"))(0.5), math.sin(0.5), 1e-12, "sin still works")
Testing.Assert_close((Compile("e(x)"))(1), math.exp(1), 1e-12, "e(x) still works")
Testing.Assert_close((Compile("e**x"))(1), math.exp(1), 1e-12, "bare e still works")
Testing.Assert_close((Compile("log_2(x)"))(8), 3, 1e-12, "log_2(x) still works")
Testing.Assert_close((Compile("sqrt(x)"))(4), 2, 1e-12, "sqrt still works")

-- ===== Bracket balancing (pre-existing behavior, locked in with a test
-- for the first time) =====

Testing.Assert_error(function() return Compile("sin((x)") end, "unclosed")
Testing.Assert_error(function() return Compile("sin(x))") end, "unmatched closing")
Testing.Assert_error(function() return Compile("sin{x)") end, "mismatched brackets")

-- ===== Compile_expression's own input validation =====

Testing.Assert_error(function() return Compile("") end, "Compile_Expression: expr")
Testing.Assert_error(function() return Compile("   ") end, "Compile_Expression: expr")
Testing.Assert_error(function() return Compile(42) end, "Compile_Expression: expr")
