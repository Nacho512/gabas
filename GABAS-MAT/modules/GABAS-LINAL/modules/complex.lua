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

-- Private tolerance for the division-by-zero guard below. Deliberately NOT
-- shared with Core's EPSILON: Core.cabs needs Is_complex from this module,
-- so if this module also required Core (just for EPSILON), the two modules
-- would require each other -- a circular dependency. Keeping this constant
-- private and duplicated here (same value, 1e-9) is a small, honest
-- tradeoff to keep the dependency graph a clean one-way arrow: Core depends
-- on Complex, never the reverse.

local DIV_EPSILON = 1e-9

local Complex_mt = {}
Complex_mt.__index = Complex_mt

local function is_complex(x)
    return type(x) == "table" and getmetatable(x) == Complex_mt
end

local function Complex(re, im)
    re, im = re or 0, im or 0
    assert(type(re) == "number" and type(im) == "number", "Complex: re and im must be numbers.")
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
    local denom = b.re * b.re + b.im * b.im
    assert(denom > DIV_EPSILON, "Complex: division by zero.")
    return Complex((a.re * b.re + a.im * b.im) / denom, (a.im * b.re - a.re * b.im) / denom)
end

function Complex_mt.__unm(a)
    return Complex(-a.re, -a.im)
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
