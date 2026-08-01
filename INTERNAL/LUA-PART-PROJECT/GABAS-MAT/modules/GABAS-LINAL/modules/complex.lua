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

-- Operator metamethods (__add.._tostring) and instance methods
-- (conj/abs/arg) for Complex_mt -- these are what let real arithmetic gain
-- Complex support "for free" throughout the rest of the library.
function Complex_mt.__add(a, b)
    a, b = to_complex(a), to_complex(b)
    return Complex(a.re + b.re, a.im + b.im)
end

function Complex_mt.__sub(a, b)
    a, b = to_complex(a), to_complex(b)
    return Complex(a.re - b.re, a.im - b.im)
end

function Complex_mt.__mul(a, b)
    a, b = to_complex(a), to_complex(b)
    return Complex(a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re)
end

function Complex_mt.__div(a, b)
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
    assert(r > 0, "Complex: 0 cannot be raised to a non-positive-integer power.")
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
}
