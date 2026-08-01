-- [Dual]
-- Dual numbers for forward-mode automatic differentiation of a REAL
-- variable: a Dual value is a table {value=..., derivative=...} with a
-- metatable implementing +, -, *, /, ^, unary -, ==, and tostring.
-- `value` and `derivative` may themselves be plain reals OR Complex values
-- (GABAS-LINAL's Complex already promotes a plain-number operand
-- automatically) -- that is what lets a REAL evaluation point carry a
-- COMPLEX-valued function and derivative, f : R -> C, which is exactly
-- the case Compile_expression needs (its `x` is always a real point).
--
-- Claude: this module deliberately does NOT support differentiating with
-- respect to a complex variable itself (Jacobian/Wirtinger modes) -- that
-- is a meaningfully different problem (holomorphicity, branch cuts), and
-- nothing here needs it, since `x` is always real. That is planned as its
-- own separate module for complex-variable calculus, not folded in here.
--
-- Because GABAS_CALC_ONE_VAR's dispatchers make every primitive
-- (sin/cos/exp/ln/sqrt/abs/...) recognize a Dual argument the same way
-- they already recognize a Complex one, differentiating a compiled
-- expression is just calling the SAME function produced by
-- Compile_expression with Dual(x0, 1) instead of a plain number -- the
-- derivative propagates through the existing compiled Lua bytecode via
-- these operator overloads, exactly how Complex already propagates
-- through GABAS-LINAL's matrix code "for free".
local Complex = require("GABAS-LINAL.modules.complex")

local Dual_mt = {}
Dual_mt.__index = Dual_mt

local function is_dual(x)
    return type(x) == "table" and getmetatable(x) == Dual_mt
end

local function Dual(value, derivative)
    derivative = derivative or 0
    return setmetatable({value = value, derivative = derivative}, Dual_mt)
end

-- A plain constant (real or Complex) has zero derivative -- differentiating
-- a number that never depends on x always gives 0.
local function to_dual(x)
    if is_dual(x) then return x end
    return Dual(x, 0)
end

-- Claude: `u`/`v` below may be plain reals or Complex -- +, -, *, / on them
-- already dispatch correctly either way through Complex's own operator
-- overloads (a plain number stays plain arithmetic; a Complex operand
-- promotes the other side automatically), so these rules need no type
-- check of their own.
function Dual_mt.__add(a, b)
    a, b = to_dual(a), to_dual(b)
    return Dual(a.value + b.value, a.derivative + b.derivative)
end

function Dual_mt.__sub(a, b)
    a, b = to_dual(a), to_dual(b)
    return Dual(a.value - b.value, a.derivative - b.derivative)
end

-- Product rule: (uv)' = u'v + uv'.
function Dual_mt.__mul(a, b)
    a, b = to_dual(a), to_dual(b)
    local u, du = a.value, a.derivative
    local v, dv = b.value, b.derivative
    return Dual(u * v, du * v + u * dv)
end

-- Quotient rule: (u/v)' = (u'v - uv') / v^2.
function Dual_mt.__div(a, b)
    a, b = to_dual(a), to_dual(b)
    local u, du = a.value, a.derivative
    local v, dv = b.value, b.derivative
    assert(not (type(v) == "number" and v == 0) and not (Complex.Is_complex(v) and v.re == 0 and v.im == 0),
        "Dual: division by zero while propagating a derivative.")
    return Dual(u / v, (du * v - u * dv) / (v * v))
end

function Dual_mt.__unm(a)
    return Dual(-a.value, -a.derivative)
end

-- Claude: general_pow(u, v) computes u^v for the VALUE slot only (no
-- derivative), used both directly (integer-exponent branch) and inside
-- the general power rule below. Negative or Complex bases are promoted to
-- Complex (same policy as Calc_one_var_sqrt/ln's auto-promotion) so this
-- never silently returns nan the way a bare native `^` would.
local function general_pow(u, v)
    if Complex.Is_complex(u) or Complex.Is_complex(v) or (type(u) == "number" and u < 0) then
        local cu = Complex.Is_complex(u) and u or Complex.Complex(u, 0)
        return cu ^ v
    end
    return u ^ v
end

local function general_ln(u)
    if Complex.Is_complex(u) then
        return Complex.Ln(u)
    end
    if u < 0 then
        return Complex.Ln(Complex.Complex(u, 0))
    end
    assert(u ~= 0, "Dual: ln(0) is undefined (encountered while differentiating a general power).")
    return math.log(u)
end

-- Claude: two branches, matching the document's guidance that integer
-- powers "should be handled directly rather than unnecessarily rewriting
-- [them] using logarithms" (more stable, and well-defined for ANY real
-- base, including negative or zero, unlike the log-based general rule).
--
-- Integer exponent, constant (not itself Dual) -- by far the most common
-- case (x^7, u^3, ...): (u^n)' = n*u^(n-1)*u'. n=0 is special-cased
-- (u^0 = 1 identically, so its derivative is exactly 0 everywhere, even
-- where u = 0 -- the direct formula would otherwise compute 0 * u^(-1),
-- which is 0 * inf = nan right where the answer is actually the simplest
-- possible case).
--
-- General case (a non-integer real exponent, or an exponent that is
-- itself Dual, i.e. depends on x too, e.g. "x**x"): w = u^v,
-- w' = w*(v'*ln(u) + v*u'/u), the general power rule. Singular at u = 0,
-- so that is rejected explicitly rather than left to propagate a nan.
function Dual_mt.__pow(a, b)
    a = to_dual(a)
    local u, du = a.value, a.derivative
    if type(b) == "number" and b == math.floor(b) and not is_dual(b) then
        local n = b
        if n == 0 then
            return Dual(1, 0)
        end
        -- Native `^` here, deliberately NOT general_pow: an INTEGER power
        -- of a negative real base is perfectly well-defined in reals
        -- ((-2)^3 = -8, confirmed against lua5.4), unlike a non-integer
        -- one -- general_pow's auto-promotion-to-Complex policy is only
        -- correct for the general branch below, not this one. `u^n` still
        -- dispatches to Complex's own __pow automatically when u already
        -- is Complex, via the same operator.
        return Dual(u ^ n, n * (u ^ (n - 1)) * du)
    end
    b = to_dual(b)
    local v, dv = b.value, b.derivative
    assert(not (type(u) == "number" and u == 0) and not (Complex.Is_complex(u) and u.re == 0 and u.im == 0),
        "Dual: the general power rule is singular at a base of 0 (only a constant integer exponent is well-defined there).")
    local w = general_pow(u, v)
    local w_deriv = w * (dv * general_ln(u) + v * (du / u))
    return Dual(w, w_deriv)
end

function Dual_mt.__eq(a, b)
    return a.value == b.value and a.derivative == b.derivative
end

function Dual_mt.__tostring(a)
    return "(" .. tostring(a.value) .. ", " .. tostring(a.derivative) .. ")"
end

return {
    Dual = Dual,
    Is_dual = is_dual,
}
