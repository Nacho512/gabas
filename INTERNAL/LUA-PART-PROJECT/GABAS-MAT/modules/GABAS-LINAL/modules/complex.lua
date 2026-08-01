-- [Complex]
-- Complex numbers: a Complex value is a table {re=..., im=...} with a
-- metatable implementing +, -, *, /, unary -, ==, and tostring. Because
-- Matrix/Tensor/Vector/Determinant/Systems/Eigen/FFT are all written using
-- plain +, -, *, / on their entries (no math.* calls baked in), they gain
-- complex-number support "for free" through these operators -- Lua looks up
-- __add/__sub/__mul/__div on whichever operand is a table, and each
-- metamethod promotes a plain number operand to Complex(x, 0) automatically.
-- A matrix/vector with only plain-number entries is completely unaffected.
--
-- This module is intentionally self-contained (no requires of its own):
-- Complex numbers are the most foundational concept in the whole library,
-- so nothing here should ever need to depend on Core or any domain module --
-- if it did, that dependency would almost certainly be pointing backwards.

-- Codex: division rejects only an exact zero denominator. Approximate-zero
-- policy belongs to callers because a fixed threshold would incorrectly
-- reject legitimate, representable small complex values.

local Complex_mt = {}
Complex_mt.__index = Complex_mt

local function is_complex(x)
    return type(x) == "table" and getmetatable(x) == Complex_mt
end

local function Complex(re, im)
    re, im = re or 0, im or 0
    assert(type(re) == "number" and type(im) == "number" and
        re == re and im == im and re ~= math.huge and re ~= -math.huge and
        im ~= math.huge and im ~= -math.huge,
        "Complex: re and im must be finite numbers.")
    return setmetatable({re = re, im = im}, Complex_mt)
end

local function to_complex(x)
    if is_complex(x) then return x end
    assert(type(x) == "number", "Complex: expected a number or a Complex value.")
    return Complex(x, 0)
end

-- Claude: Lua only ever dispatches a binary operator to ONE operand's
-- metamethod (the left operand's, if it has one, else the right's) --
-- so "Complex_value * Dual_value" always lands in Complex_mt.__mul, never
-- Dual's, even though Dual (GABAS-CALC-ONE-VAR's forward-mode automatic
-- differentiation number) is the one that actually knows how to combine
-- the two. Without this, that combination would fail outright inside
-- to_complex, which has no idea what a Dual is -- correctly so: Complex
-- must stay self-contained and never require a module built on top of it.
--
-- The fix is to detect a genuinely FOREIGN operand (a table that isn't
-- Complex but defines this same metamethod) and hand the whole operation
-- to IT instead -- Complex never needs to know Dual (or anything else)
-- exists, it just declines to handle a type it doesn't recognize when
-- something else clearly can. Dual's own arithmetic then treats the
-- Complex operand as a plain constant (zero derivative), which is exactly
-- right.
local function foreign_metamethod(x, op)
    if type(x) == "table" and not is_complex(x) then
        local mt = getmetatable(x)
        if mt then return mt[op] end
    end
    return nil
end

-- Operator metamethods (__add.._tostring) and instance methods
-- (conj/abs/arg) for Complex_mt -- these are what let real arithmetic gain
-- Complex support "for free" throughout the rest of the library.
function Complex_mt.__add(a, b)
    local foreign = foreign_metamethod(a, "__add") or foreign_metamethod(b, "__add")
    if foreign then return foreign(a, b) end
    a, b = to_complex(a), to_complex(b)
    return Complex(a.re + b.re, a.im + b.im)
end

function Complex_mt.__sub(a, b)
    local foreign = foreign_metamethod(a, "__sub") or foreign_metamethod(b, "__sub")
    if foreign then return foreign(a, b) end
    a, b = to_complex(a), to_complex(b)
    return Complex(a.re - b.re, a.im - b.im)
end

function Complex_mt.__mul(a, b)
    local foreign = foreign_metamethod(a, "__mul") or foreign_metamethod(b, "__mul")
    if foreign then return foreign(a, b) end
    a, b = to_complex(a), to_complex(b)
    return Complex(a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re)
end

function Complex_mt.__div(a, b)
    local foreign = foreign_metamethod(a, "__div") or foreign_metamethod(b, "__div")
    if foreign then return foreign(a, b) end
    a, b = to_complex(a), to_complex(b)
    assert(b.re ~= 0 or b.im ~= 0, "Complex: division by zero.")
    if math.abs(b.re) >= math.abs(b.im) then
        local ratio = b.im / b.re
        local denom = b.re + b.im * ratio
        return Complex((a.re + a.im * ratio) / denom, (a.im - a.re * ratio) / denom)
    end
    local ratio = b.re / b.im
    local denom = b.im + b.re * ratio
    return Complex((a.re * ratio + a.im) / denom, (a.im * ratio - a.re) / denom)
end

function Complex_mt.__unm(a)
    return Complex(-a.re, -a.im)
end

-- Claude: a real integer exponent goes through exact repeated squaring
-- (so e.g. j^2 comes out to exactly -1, not cos/sin round-off) -- every
-- other case (a non-integer real exponent, or a genuinely Complex one,
-- e.g. e^j or j^j) falls through to the general definition of complex
-- exponentiation, z^w = exp(w * Log(z)), using the principal complex
-- logarithm Log(z) = ln|z| + i*arg(z). This is the only way a Complex
-- exponent is defined at all, and it also subsumes the old real,
-- non-integer case (De Moivre's formula is exactly what this reduces to
-- when w is real), so there is no longer a separate branch for it.
function Complex_mt.__pow(a, b)
    local foreign = foreign_metamethod(a, "__pow") or foreign_metamethod(b, "__pow")
    if foreign then return foreign(a, b) end
    a = to_complex(a)
    if type(b) == "number" and b == math.floor(b) then
        local n = math.abs(b)
        local result = Complex(1, 0)
        local base = a
        while n > 0 do
            if n % 2 == 1 then result = result * base end
            base = base * base
            n = math.floor(n / 2)
        end
        if b < 0 then
            return to_complex(1) / result
        end
        return result
    end
    b = to_complex(b)
    local r = a:abs()
    if r == 0 then
        -- Claude: 0^v is 0 by the standard real-analysis convention
        -- whenever the exponent's real part is positive (continuous with
        -- lim r->0+ of r^v) -- e.g. Sqrt(0) = 0^0.5 must be 0, not an
        -- error (this surfaced via Asin(1), whose formula computes
        -- Sqrt(1-1^2)). Only a nonpositive real part is genuinely
        -- undefined (0^0 itself is handled separately, by the exact
        -- integer branch above, before this general branch is ever
        -- reached).
        assert(b.re > 0, "Complex: 0 cannot be raised to a power with a nonpositive real part.")
        return Complex(0, 0)
    end
    local ln_r = math.log(r)
    local theta = a:arg()
    local re_exp = b.re * ln_r - b.im * theta
    local im_exp = b.re * theta + b.im * ln_r
    local mag = math.exp(re_exp)
    return Complex(mag * math.cos(im_exp), mag * math.sin(im_exp))
end

function Complex_mt.__eq(a, b)
    return a.re == b.re and a.im == b.im
end

function Complex_mt.__tostring(a)
    if a.im == 0 then
        return tostring(a.re)
    end
    local sign = a.im < 0 and "-" or "+"
    return tostring(a.re) .. sign .. tostring(math.abs(a.im)) .. "i"
end

function Complex_mt.conj(a)
    return Complex(a.re, -a.im)
end

function Complex_mt.abs(a)
    return math.sqrt(a.re * a.re + a.im * a.im)
end

function Complex_mt.arg(a)
    return math.atan(a.im, a.re)
end

-- Claude: private, real-only sinh/cosh -- not exported. Needed just as
-- building blocks for Sin/Cos/Sinh/Cosh below (their closed forms mix the
-- real trig and real hyperbolic functions of z's components). Built from
-- the exponential definition rather than relying on math.sinh/cosh, same
-- reasoning as Calc_one_var_sinh/cosh in GABAS-CALC-ONE-VAR: those were
-- made optional in Lua 5.3+ and this module cannot depend on that module
-- (or vice versa) without pointing the dependency backwards.
local function real_sinh(t) return (math.exp(t) - math.exp(-t)) / 2 end
local function real_cosh(t) return (math.exp(t) + math.exp(-t)) / 2 end

-- Exp/Ln/Sqrt/Sin/Cos/Tan/Sinh/Cosh/Tanh: the elementary functions,
-- generalized to a Complex argument (a plain real number is accepted too,
-- coerced via to_complex, so these work standalone as well as from a
-- dispatcher that only calls them once a value is already known Complex).
-- Ln and Sqrt use the SAME principal branch as __pow above (Log(z) =
-- ln|z| + i*arg(z), with arg via atan2-style math.atan(im, re)) -- Sqrt is
-- literally just the z^0.5 case of __pow, reused rather than re-derived.
local function Exp(z)
    z = to_complex(z)
    local mag = math.exp(z.re)
    return Complex(mag * math.cos(z.im), mag * math.sin(z.im))
end

local function Ln(z)
    z = to_complex(z)
    local r = z:abs()
    assert(r > 0, "Complex: ln(0) is undefined.")
    return Complex(math.log(r), z:arg())
end

local function Sqrt(z)
    return to_complex(z) ^ 0.5
end

local function Sin(z)
    z = to_complex(z)
    return Complex(math.sin(z.re) * real_cosh(z.im), math.cos(z.re) * real_sinh(z.im))
end

local function Cos(z)
    z = to_complex(z)
    return Complex(math.cos(z.re) * real_cosh(z.im), -math.sin(z.re) * real_sinh(z.im))
end

local function Tan(z)
    return Sin(z) / Cos(z)
end

local function Sinh(z)
    z = to_complex(z)
    return Complex(real_sinh(z.re) * math.cos(z.im), real_cosh(z.re) * math.sin(z.im))
end

local function Cosh(z)
    z = to_complex(z)
    return Complex(real_cosh(z.re) * math.cos(z.im), real_sinh(z.re) * math.sin(z.im))
end

local function Tanh(z)
    return Sinh(z) / Cosh(z)
end

-- Claude: Asin/Acos/Atan via the standard logarithmic definitions --
-- asin(z) = -i*ln(iz + sqrt(1-z^2)), acos(z) = -i*ln(z + i*sqrt(1-z^2)),
-- atan(z) = (i/2)*ln((1-iz)/(1+iz)) -- reusing Ln and Sqrt (and their
-- already-principal-branch behavior) rather than re-deriving branch
-- handling from scratch. Singular where the formula's own denominator or
-- logarithm argument vanishes (Atan at z = +-i in particular, matching
-- the well-known branch points of arctangent); not specially guarded
-- here, same as Ln(0) or Sqrt's own domain limits -- the underlying
-- arithmetic already fails loudly (division by zero, ln of zero) rather
-- than silently.
local IMAGINARY_UNIT = Complex(0, 1)
local ONE = Complex(1, 0)

local function Asin(z)
    z = to_complex(z)
    return -IMAGINARY_UNIT * Ln(IMAGINARY_UNIT * z + Sqrt(ONE - z * z))
end

local function Acos(z)
    z = to_complex(z)
    return -IMAGINARY_UNIT * Ln(z + IMAGINARY_UNIT * Sqrt(ONE - z * z))
end

local function Atan(z)
    z = to_complex(z)
    return (IMAGINARY_UNIT / 2) * Ln((ONE - IMAGINARY_UNIT * z) / (ONE + IMAGINARY_UNIT * z))
end

-- Elementwise conjugate. Accepts either a vector (flat table) or a
-- matrix (table of row-tables); plain-number entries pass through
-- unchanged.
local function Conjugate(x)
    local C = {}
    if type(x[1]) == "table" and not is_complex(x[1]) then
        for i = 1, #x do
            local row, row_c = x[i], {}
            for j = 1, #row do
                row_c[j] = is_complex(row[j]) and row[j]:conj() or row[j]
            end
            C[i] = row_c
        end
    else
        for i = 1, #x do
            C[i] = is_complex(x[i]) and x[i]:conj() or x[i]
        end
    end
    return C
end

return {
    Complex = Complex,
    Is_complex = is_complex,
    Conjugate = Conjugate,
    Exp = Exp,
    Ln = Ln,
    Sqrt = Sqrt,
    Sin = Sin,
    Cos = Cos,
    Tan = Tan,
    Sinh = Sinh,
    Cosh = Cosh,
    Tanh = Tanh,
    Asin = Asin,
    Acos = Acos,
    Atan = Atan,
}
