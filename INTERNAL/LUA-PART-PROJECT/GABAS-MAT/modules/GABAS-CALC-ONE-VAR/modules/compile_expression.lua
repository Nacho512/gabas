-- [Compile Expression]
-- Compile_Expression -- turns a string the user typed in a friendly,
-- calculator-like syntax into a real, callable Lua function of one
-- variable x. This is the foundational piece dif_numeric.lua and
-- integration_numeric.lua both build on: they operate on whatever
-- function this module produces, without needing to know anything about
-- how it was parsed.
--
-- Claude: every trig function here (sin/cos/tan and their Complex-aware
-- versions) works in RADIANS, exclusively, with no implicit conversion
-- anywhere in the pipeline -- deliberately, not just because that's what
-- math.sin/cos/tan happen to do. This is the one convention every
-- differentiation/integration identity in calculus assumes: d/dx sin(x)
-- = cos(x) is only true when x is in radians -- in degrees, the chain
-- rule silently introduces a pi/180 factor into every derivative that
-- touches a trig function, since sin(x_deg) really means
-- sin(x_deg * pi/180). This module must never blur that: internally,
-- "x" is always radians, full stop. Degree INPUT is offered separately
-- and explicitly (sind/cosd/tand below), as pure syntax sugar that
-- converts to radians before ever reaching sin/cos/tan -- so the
-- pi/180 factor stays a normal, visible consequence of the chain rule
-- once derivatives are built, rather than something hidden inside the
-- trig functions themselves.
local Complex = require("GABAS-LINAL.modules.complex")
local Dual = require("GABAS-CALC-ONE-VAR.modules.dual")
local Core = require("GABAS-CALC-ONE-VAR.modules.core")

-- Claude: whitespace is insignificant everywhere in a math expression (no
-- string literals ever appear here), so every step downstream can assume a
-- fully compact string -- stripped once, up front, rather than have every
-- pattern below account for optional gaps around parens, "**", unary "-",
-- etc. Negative exponents (e.g. "x**-2") need no special handling beyond
-- this: Lua's own grammar already accepts unary "-" as the operand of "^".
local function strip_whitespace(expr)
    return (expr:gsub("%s+", ""))
end

-- Recognized opening brackets, all meant as pure grouping -- "(", "{", and
-- "[" carry no special meaning here, unlike in Lua itself (table
-- constructors, indexing). What they must NOT do is mismatch: opening
-- with "{" and closing with "]" is rejected, even though both eventually
-- become "(" / ")" once normalized. Checked below with a stack that
-- tracks bracket TYPE, not just a running open/close count.
local OPEN = { ["("] = true, ["{"] = true, ["["] = true }
local CLOSE_FOR_OPEN = { ["("] = ")", ["{"] = "}", ["["] = "]" }

-- Claude: validates that every bracket in `expr` is balanced AND correctly
-- typed (a "{" can only be closed by "}", never by ")" or "]"), then
-- returns a new string with every bracket normalized to plain "(" / ")" --
-- so every step after this one only ever has to deal with one bracket
-- style, no matter which ones the user actually typed.
local function normalize_brackets(expr)
    local stack = {}
    local out = {}
    for i = 1, #expr do
        local c = expr:sub(i, i)
        if OPEN[c] then
            stack[#stack + 1] = c
            out[#out + 1] = "("
        elseif c == ")" or c == "}" or c == "]" then
            local opener = stack[#stack]
            assert(opener,
                "Compile_Expression: unmatched closing '" .. c .. "' at position " ..
                i .. " -- no bracket was open to close.")
            assert(CLOSE_FOR_OPEN[opener] == c,
                "Compile_Expression: mismatched brackets -- opened with '" .. opener ..
                "' but closed with '" .. c .. "' at position " .. i ..
                " (must close with '" .. CLOSE_FOR_OPEN[opener] .. "').")
            stack[#stack] = nil
            out[#out + 1] = ")"
        else
            out[#out + 1] = c
        end
    end
    assert(#stack == 0,
        "Compile_Expression: unclosed '" .. tostring(stack[#stack]) ..
        "' -- every opened bracket must be closed.")
    return table.concat(out)
end

-- log_N(...) -> a base-N logarithm, N given as a literal integer suffix
-- (log_3(...), log_10(...), etc.), matching the subscript notation from
-- Nacho's handwritten example. Rewritten as a PREFIX substitution --
-- "log_3(" becomes "Calc_one_var_log_base(3," -- rather than trying to
-- relocate the "3" to after the argument (which would need to find the
-- matching close-paren by hand). Lua's own parser, reading the result,
-- finds that matching paren for free, since after this rewrite it's just
-- an ordinary function call.
local function translate_log_base(expr)
    return (expr:gsub("log_(%d+)%(", "Calc_one_var_log_base(%1,"))
end

-- root_N(...) -> the N-th root of (...), N given as a literal (possibly
-- decimal) index suffix (root_3(x), root_6.13(x), ...), the same
-- prefix-substitution trick as log_N above: "root_6.13(" becomes
-- "Root_any(6.13,", and Lua's own parser finds the matching close-paren.
local function translate_root_index(expr)
    return (expr:gsub("root_([%d%.]+)%(", "Root_any(%1,"))
end

-- e(...) -> Euler's number raised to the argument (e.g. e(x+1) means
-- e^(x+1)), a prefix substitution exactly like log_N above -- Lua's own
-- parser finds the matching close-paren for free once this is rewritten
-- as an ordinary function call. Routed through the Calc_one_var_exp
-- dispatcher (not straight to math.exp) so e(j) and friends work the same
-- as exp(j).
local function translate_e_call(expr)
    return (expr:gsub("%f[%a]e%(", "Calc_one_var_exp("))
end

-- Claude: bare "e" as a standalone constant (so "e**{x+1}" -- after
-- bracket/** normalization, "e^(x+1)" -- means e^(x+1), not a syntax
-- error). This can't be a simple frontier-pattern gsub like the others:
-- Lua's OWN numeric literals use a trailing "e"/"E" for scientific
-- notation (5e3, 1.5e-10), and that "e" must never be touched. The one
-- reliable signal that a lone "e" is genuinely a standalone identifier
-- (ours) rather than a numeric-literal exponent marker is that a
-- literal's "e" is always immediately preceded by a digit; ours never
-- is. Scanned by hand, not gsub, so each candidate "e" can be checked
-- against its actual neighboring characters, not just a pattern class.
-- Must run AFTER translate_e_call above, so "e(" cases are already
-- consumed and never reach this more general rule.
local function translate_bare_e(expr)
    local out = {}
    local n = #expr
    for i = 1, n do
        local c = expr:sub(i, i)
        if c == "e" then
            local before = i > 1 and expr:sub(i - 1, i - 1) or ""
            local after = i < n and expr:sub(i + 1, i + 1) or ""
            if not before:match("%w") and not after:match("%w") then
                out[#out + 1] = "Calc_one_var_e"
            else
                out[#out + 1] = c
            end
        else
            out[#out + 1] = c
        end
    end
    return table.concat(out)
end

-- Claude: word-boundary substitutions (Lua's %f frontier pattern) so
-- "sin" only matches as a genuinely standalone identifier -- never as a
-- substring of "asin" or "sinh". This is what lets every entry below be
-- applied in any order without one clobbering another.
--
-- Every entry routes to a dispatcher (Calc_one_var_*, defined below)
-- rather than straight to math.* -- each one accepts a plain real number
-- (unchanged behavior, still returns a plain number), a Complex value
-- (once "j" is in play), or a Dual value (once Derivative_at_point is
-- propagating a derivative through it), dispatching to the matching
-- formula in each case.
--
-- sind/cosd/tand are the degree-input convenience versions -- see the note
-- above local Calc_one_var_sind for why this is safe to bolt on as pure
-- syntax sugar without touching what sin/cos/tan themselves mean.
local FUNCTION_ALIASES = {
    sin = "Calc_one_var_sin", cos = "Calc_one_var_cos", tan = "Calc_one_var_tan",
    sind = "Calc_one_var_sind", cosd = "Calc_one_var_cosd", tand = "Calc_one_var_tand",
    asin = "Calc_one_var_asin", acos = "Calc_one_var_acos", atan = "Calc_one_var_atan",
    sinh = "Calc_one_var_sinh", cosh = "Calc_one_var_cosh", tanh = "Calc_one_var_tanh",
    abs = "Calc_one_var_abs", sqrt = "Calc_one_var_sqrt", exp = "Calc_one_var_exp", ln = "Calc_one_var_ln",
}

local function translate_function_names(expr)
    for name, replacement in pairs(FUNCTION_ALIASES) do
        expr = expr:gsub("%f[%a]" .. name .. "%f[%A]", replacement)
    end
    return expr
end

-- Claude: "j" as the imaginary unit (engineering convention, sqrt(-1)),
-- e.g. "3+4*j". Unlike bare "e", this needs no careful hand-scan: "j" is
-- not part of Lua's own numeric-literal grammar (no exponent marker, no
-- hex-float suffix uses it), so a plain frontier-pattern gsub is safe.
local function translate_imaginary_unit(expr)
    return (expr:gsub("%f[%a]j%f[%A]", "Calc_one_var_j"))
end

-- Claude: sinh/cosh/tanh from their exponential definitions, NOT Lua's
-- own math.sinh/cosh/tanh -- those were made optional (LUA_COMPAT_5_2) as
-- of Lua 5.3 and aren't guaranteed to exist in every build, including
-- whatever Lua LuaTeX embeds, which is the actual target for this module.
-- Defining these ourselves means this never depends on a compile flag we
-- don't control. Private (real-only): the public sinh/cosh/tanh below are
-- the dispatchers that call these for a real argument.
local function real_sinh(x) return (math.exp(x) - math.exp(-x)) / 2 end
local function real_cosh(x) return (math.exp(x) + math.exp(-x)) / 2 end
local function real_tanh(x) return real_sinh(x) / real_cosh(x) end

-- Claude: real_fn is used for a plain number (unchanged behavior, still
-- returns a plain number); complex_fn (from GABAS-LINAL's Complex module)
-- is used once an argument is actually Complex, i.e. some upstream "j"
-- was involved; dual_fn is used once an argument is actually Dual, i.e.
-- Derivative_at_point is propagating a derivative through this call. This
-- is the shared shape for every dispatcher below.
local function dispatch3(real_fn, complex_fn, dual_fn)
    return function(z)
        if Dual.Is_dual(z) then
            return dual_fn(z)
        end
        if Complex.Is_complex(z) then
            return complex_fn(z)
        end
        return real_fn(z)
    end
end

-- Claude: forward-declared so the dual_fn rules just below can close over
-- them by upvalue (e.g. sin's rule needs cos) even though they're only
-- ASSIGNED after every rule is defined -- by the time any of these is
-- actually called, every upvalue already has its final function.
local Calc_one_var_sin, Calc_one_var_cos, Calc_one_var_tan
local Calc_one_var_sinh, Calc_one_var_cosh, Calc_one_var_tanh

-- Claude: each dual_fn is the chain rule, g(u)' = g'(u)*u', reusing the
-- ALREADY-DISPATCHED real/Complex primitives above for both g(u) and the
-- derivative-rule term g'(u) -- u itself is never Dual here (it's the
-- .value of an outer Dual), so calling e.g. Calc_one_var_cos(u) safely
-- lands in the real or Complex branch, never recurses back into a dual_fn.
local function sin_dual(z)
    local u, du = z.value, z.derivative
    return Dual.Dual(Calc_one_var_sin(u), Calc_one_var_cos(u) * du)
end
local function cos_dual(z)
    local u, du = z.value, z.derivative
    return Dual.Dual(Calc_one_var_cos(u), -(Calc_one_var_sin(u) * du))
end
local function tan_dual(z)
    local u, du = z.value, z.derivative
    local c = Calc_one_var_cos(u)
    return Dual.Dual(Calc_one_var_tan(u), du / (c * c))
end
local function sinh_dual(z)
    local u, du = z.value, z.derivative
    return Dual.Dual(Calc_one_var_sinh(u), Calc_one_var_cosh(u) * du)
end
local function cosh_dual(z)
    local u, du = z.value, z.derivative
    return Dual.Dual(Calc_one_var_cosh(u), Calc_one_var_sinh(u) * du)
end
local function tanh_dual(z)
    local u, du = z.value, z.derivative
    local c = Calc_one_var_cosh(u)
    return Dual.Dual(Calc_one_var_tanh(u), du / (c * c))
end

Calc_one_var_sin = dispatch3(math.sin, Complex.Sin, sin_dual)
Calc_one_var_cos = dispatch3(math.cos, Complex.Cos, cos_dual)
Calc_one_var_tan = dispatch3(math.tan, Complex.Tan, tan_dual)
Calc_one_var_sinh = dispatch3(real_sinh, Complex.Sinh, sinh_dual)
Calc_one_var_cosh = dispatch3(real_cosh, Complex.Cosh, cosh_dual)
Calc_one_var_tanh = dispatch3(real_tanh, Complex.Tanh, tanh_dual)

-- Claude: sind/cosd/tand -- degree-input convenience, but the ONLY thing
-- they do is multiply by DEG_TO_RAD before handing off to the real
-- (radians-only) dispatchers above; nothing about sin/cos/tan themselves
-- changes. z*DEG_TO_RAD works unchanged whether z is a plain real, a
-- Complex, or a Dual value (each operator already promotes a plain-number
-- operand), so these need no type check of their own -- and because the
-- conversion is ordinary multiplication rather than something
-- special-cased, differentiating through it picks up the pi/180 factor
-- for free via the ordinary chain rule, exactly like differentiating any
-- other composition.
local DEG_TO_RAD = math.pi / 180
local function Calc_one_var_sind(z) return Calc_one_var_sin(z * DEG_TO_RAD) end
local function Calc_one_var_cosd(z) return Calc_one_var_cos(z * DEG_TO_RAD) end
local function Calc_one_var_tand(z) return Calc_one_var_tan(z * DEG_TO_RAD) end

-- Claude: exp/ln/sqrt need one more branch than a plain real/Complex
-- dispatch: a NEGATIVE real argument to ln or sqrt has no real result
-- (today that would silently come back as nan from math.log/math.sqrt),
-- but it has a perfectly well-defined Complex one (sqrt(-4) = 2i,
-- ln(-1) = i*pi) now that Complex.Sqrt/Ln exist -- so a negative real is
-- auto-promoted to Complex rather than left to fail silently. Forward-
-- declared for the same self-reference reason as sin/cos/... above (each
-- dual_fn calls the real/Complex-dispatched version of itself, on u).
local Calc_one_var_exp, Calc_one_var_ln, Calc_one_var_sqrt, Calc_one_var_abs

-- d/du[e^u] = e^u -- reuses new_value (already e^u) rather than
-- recomputing it.
local function exp_dual(z)
    local u, du = z.value, z.derivative
    local new_value = Calc_one_var_exp(u)
    return Dual.Dual(new_value, new_value * du)
end

Calc_one_var_exp = function(z)
    if Dual.Is_dual(z) then return exp_dual(z) end
    if Complex.Is_complex(z) then return Complex.Exp(z) end
    return math.exp(z)
end

-- d/du[ln(u)] = u'/u -- valid unchanged whether u ended up real or
-- Complex-promoted: promoting to the principal branch only adds a
-- CONSTANT i*pi, which vanishes under differentiation, so the formula
-- itself never needs to branch on which case produced the value.
local function ln_dual(z)
    local u, du = z.value, z.derivative
    return Dual.Dual(Calc_one_var_ln(u), du / u)
end

Calc_one_var_ln = function(z)
    if Dual.Is_dual(z) then return ln_dual(z) end
    if Complex.Is_complex(z) then return Complex.Ln(z) end
    if z < 0 then return Complex.Ln(Complex.Complex(z, 0)) end
    return math.log(z)
end

-- d/du[sqrt(u)] = 1/(2*sqrt(u)) -- reuses new_value (already sqrt(u)),
-- same "the branch is just a constant offset" reasoning as ln above.
local function sqrt_dual(z)
    local u, du = z.value, z.derivative
    local new_value = Calc_one_var_sqrt(u)
    return Dual.Dual(new_value, du / (2 * new_value))
end

Calc_one_var_sqrt = function(z)
    if Dual.Is_dual(z) then return sqrt_dual(z) end
    if Complex.Is_complex(z) then return Complex.Sqrt(z) end
    if z < 0 then return Complex.Sqrt(Complex.Complex(z, 0)) end
    return math.sqrt(z)
end

-- Claude: |x| is genuinely NOT differentiable at x = 0 (the left and
-- right derivatives disagree, -1 vs +1) -- rejected explicitly rather
-- than arbitrarily returning 0 the way a careless branch might (the
-- specification's own warning about the abs primitive). For a
-- Complex-valued u(x) -- which is still a function of a single REAL x, so
-- this stays well-defined -- |u(x)| = sqrt(re^2+im^2) is real-valued, and
-- its derivative is Re(u' * conj(u)) / |u|, away from u = 0.
local function abs_dual(z)
    local u, du = z.value, z.derivative
    if Complex.Is_complex(u) then
        local new_value = u:abs()
        assert(new_value ~= 0, "Calc_one_var_abs: not differentiable at 0.")
        local product = du * u:conj()
        local re_part = Complex.Is_complex(product) and product.re or product
        return Dual.Dual(new_value, re_part / new_value)
    end
    assert(u ~= 0, "Calc_one_var_abs: not differentiable at 0 (left and right derivatives disagree).")
    local sign = u > 0 and 1 or -1
    return Dual.Dual(math.abs(u), sign * du)
end

Calc_one_var_abs = function(z)
    if Dual.Is_dual(z) then return abs_dual(z) end
    if Complex.Is_complex(z) then return z:abs() end
    return math.abs(z)
end

-- Claude: asin/acos need the same real-domain-failure promotion as ln/
-- sqrt above -- |x|>1 has no real arcsine/arccosine, but a perfectly
-- well-defined Complex one. atan has no such restriction (every real
-- number has a real arctangent). Derivative rules reuse the ALREADY-
-- dispatched sqrt (so the same negative-argument-promotes-to-Complex
-- policy carries through automatically): d/du[asin(u)] = u'/sqrt(1-u^2),
-- d/du[acos(u)] = -u'/sqrt(1-u^2), d/du[atan(u)] = u'/(1+u^2).
local Calc_one_var_asin, Calc_one_var_acos, Calc_one_var_atan

local function asin_dual(z)
    local u, du = z.value, z.derivative
    local denom = Calc_one_var_sqrt(1 - u * u)
    return Dual.Dual(Calc_one_var_asin(u), du / denom)
end

Calc_one_var_asin = function(z)
    if Dual.Is_dual(z) then return asin_dual(z) end
    if Complex.Is_complex(z) then return Complex.Asin(z) end
    if z < -1 or z > 1 then return Complex.Asin(Complex.Complex(z, 0)) end
    return math.asin(z)
end

local function acos_dual(z)
    local u, du = z.value, z.derivative
    local denom = Calc_one_var_sqrt(1 - u * u)
    return Dual.Dual(Calc_one_var_acos(u), -(du / denom))
end

Calc_one_var_acos = function(z)
    if Dual.Is_dual(z) then return acos_dual(z) end
    if Complex.Is_complex(z) then return Complex.Acos(z) end
    if z < -1 or z > 1 then return Complex.Acos(Complex.Complex(z, 0)) end
    return math.acos(z)
end

local function atan_dual(z)
    local u, du = z.value, z.derivative
    return Dual.Dual(Calc_one_var_atan(u), du / (1 + u * u))
end

Calc_one_var_atan = function(z)
    if Dual.Is_dual(z) then return atan_dual(z) end
    if Complex.Is_complex(z) then return Complex.Atan(z) end
    return math.atan(z)
end

-- Claude: same negative-real-promotes-to-Complex policy as ln above,
-- generalized to log_base(x) = Ln(x) / Ln(base) (base is always a real,
-- positive literal from "log_N(...)" -- see translate_log_base, so base
-- itself is never Dual). Derivative: d/du[log_base(u)] = u'/(u*ln(base)),
-- same "branch is just a constant offset" reasoning as ln.
local function Calc_one_var_log_base(base, x)
    if Dual.Is_dual(x) then
        local u, du = x.value, x.derivative
        local new_value = Calc_one_var_log_base(base, u)
        return Dual.Dual(new_value, du / (u * math.log(base)))
    end
    if Complex.Is_complex(x) then
        return Complex.Ln(x) / math.log(base)
    end
    if x < 0 then
        return Complex.Ln(Complex.Complex(x, 0)) / math.log(base)
    end
    return math.log(x, base)
end

-- The N-th (possibly non-integer) root of radicand, real-valued whenever
-- one exists: directly as radicand^(1/index) when radicand >= 0, or, for a
-- negative radicand, only when index is an odd integer (the one case a
-- real root of a negative number is still well-defined, e.g. the cube
-- root of -8 is -2) -- computed as -((-radicand)^(1/index)) to avoid
-- raising a negative base to a fractional power, which Lua's `^` cannot
-- do (it would silently produce nan).
--
-- Claude: the Dual branch mirrors this SAME real-only convention (it does
-- NOT delegate to the general power rule in dual.lua, which would
-- auto-promote a negative radicand to Complex -- a different, inconsistent
-- answer from Root_any's own -8^(1/3) = -2 policy). index is always a
-- plain constant here (never Dual), since "root_N(...)" only ever
-- captures a literal number for N. A Complex radicand is not supported,
-- consistent with the non-Dual branch just below rejecting one too.
local function Root_any(index, radicand)
    Core.assert_finite_number(index, "Root_any: index")
    assert(index ~= 0, "Root_any: index cannot be zero.")
    if Dual.Is_dual(radicand) then
        local r, dr = radicand.value, radicand.derivative
        Core.assert_finite_number(r, "Root_any: radicand")
        if r > 0 then
            local new_value = r ^ (1 / index)
            return Dual.Dual(new_value, (1 / index) * (r ^ (1 / index - 1)) * dr)
        end
        assert(index == math.floor(index) and (index % 2) == 1,
            "Root_any: the derivative of the root of a negative radicand is only real-valued when index is an odd integer.")
        assert(r ~= 0, "Root_any: not differentiable at radicand = 0.")
        local new_value = -((-r) ^ (1 / index))
        return Dual.Dual(new_value, (1 / index) * ((-r) ^ (1 / index - 1)) * dr)
    end
    Core.assert_finite_number(radicand, "Root_any: radicand")
    if radicand >= 0 then
        return radicand ^ (1 / index)
    end
    assert(index == math.floor(index) and (index % 2) == 1,
        "Root_any: the root of a negative radicand is only real-valued when index is an odd integer.")
    return -((-radicand) ^ (1 / index))
end

-- Precomputed once, like math.pi, rather than recomputed on every
-- reference -- Euler's number, for the bare-"e" case (e**{x+1}).
local Calc_one_var_e = math.exp(1)

-- Precomputed once -- the imaginary unit, for the bare-"j" case (3+4*j).
local Calc_one_var_j = Complex.Complex(0, 1)

-- The environment the compiled expression actually runs in -- deliberately
-- NOT Lua's real global table (_G), and deliberately NOT the whole `math`
-- table either (see ALLOWED_IDENTIFIERS below for why). Compiled user
-- expressions only ever see the handful of helpers below; they can't
-- read or write anything else in the process, and nothing in the process
-- can be shadowed by whatever a user happens to type.
local COMPILE_ENV = {
    Calc_one_var_log_base = Calc_one_var_log_base,
    Calc_one_var_sin = Calc_one_var_sin,
    Calc_one_var_cos = Calc_one_var_cos,
    Calc_one_var_tan = Calc_one_var_tan,
    Calc_one_var_sind = Calc_one_var_sind,
    Calc_one_var_cosd = Calc_one_var_cosd,
    Calc_one_var_tand = Calc_one_var_tand,
    Calc_one_var_asin = Calc_one_var_asin,
    Calc_one_var_acos = Calc_one_var_acos,
    Calc_one_var_atan = Calc_one_var_atan,
    Calc_one_var_sinh = Calc_one_var_sinh,
    Calc_one_var_cosh = Calc_one_var_cosh,
    Calc_one_var_tanh = Calc_one_var_tanh,
    Calc_one_var_exp = Calc_one_var_exp,
    Calc_one_var_ln = Calc_one_var_ln,
    Calc_one_var_sqrt = Calc_one_var_sqrt,
    Calc_one_var_abs = Calc_one_var_abs,
    Calc_one_var_e = Calc_one_var_e,
    Calc_one_var_j = Calc_one_var_j,
    Root_any = Root_any,
}

-- Claude: the exhaustive whitelist of every bare identifier this package
-- permits in the FULLY NORMALIZED source about to be compiled (technical
-- document, section 17.1's "strict whitelist" architecture) -- used as a
-- Lua table the same way a Python dictionary would be, purely for O(1)
-- membership testing (name -> true). "x" (the one real independent
-- variable) is added by hand, since it is a function PARAMETER in the
-- compiled source, not a COMPILE_ENV key; every other allowed name is
-- read directly off COMPILE_ENV's own keys, rather than duplicated here
-- by hand -- so this can never silently drift out of sync with what the
-- compiled code can actually reach. Adding a new dispatcher to
-- COMPILE_ENV automatically whitelists it; nothing can be reachable from
-- compiled code without ALSO being whitelisted here, and nothing can be
-- whitelisted here without ALSO actually being reachable.
local ALLOWED_IDENTIFIERS = { x = true }
for key in pairs(COMPILE_ENV) do
    ALLOWED_IDENTIFIERS[key] = true
end

-- Claude: strips every Lua NUMERIC LITERAL (including a scientific-
-- notation exponent, e.g. "1e10", "6.13", ".5e-3") out of a COPY of the
-- source before scanning it for identifiers -- otherwise the exponent
-- marker in a numeral like "1e10" would be misread as a stray identifier
-- "e10", exactly the ambiguity translate_bare_e above already has to
-- navigate for the SAME reason. Tried most-specific-first (an exponent
-- form before the bare form it is a superset of), so a numeral is always
-- consumed whole, in one pass, never leaving a fragment like a lone "e10"
-- behind for the plainer pattern to mis-split.
local function mask_numeric_literals(expr)
    expr = expr:gsub("%d+%.?%d*[eE][%+%-]?%d+", "")
    expr = expr:gsub("%.%d+[eE][%+%-]?%d+", "")
    expr = expr:gsub("%d+%.?%d*", "")
    expr = expr:gsub("%.%d+", "")
    return expr
end

-- Claude: rejects a normalized expression that references any bare
-- identifier OUTSIDE ALLOWED_IDENTIFIERS -- a typo ("sni(x)"), a
-- disallowed global ("math.random()", now that "math" itself was
-- deliberately removed from COMPILE_ENV so it could never be reached
-- this way), or a Lua keyword smuggled in as if it were part of the
-- expression ("x and 1" -- "and" is lexically indistinguishable from an
-- identifier at this level, and simply is not on the whitelist). This
-- runs BEFORE load() ever sees the source, so a disallowed identifier
-- stops compilation immediately, with the exact offending name named in
-- the error, rather than being deferred to a confusing runtime failure
-- (an undefined global reads as nil, so calling it raises "attempt to
-- call a nil value") the first time some LATER, unrelated caller happens
-- to invoke the compiled function.
local function assert_allowed_identifiers(normalized, original)
    for name in mask_numeric_literals(normalized):gmatch("%a[%w_]*") do
        assert(ALLOWED_IDENTIFIERS[name],
            "Compile_Expression: \"" .. original .. "\" uses \"" .. name ..
            "\", which is not part of this package's approved expression syntax " ..
            "(not a recognized function, constant, or the variable x).")
    end
end

-- Compile_Expression(expr): turns a string like "sin(x^2) + 3" into a
-- real Lua function f(x). `**` is accepted as an alternative to Lua's
-- native `^` (both work; `**` is translated, `^` is left untouched since
-- Lua already understands it). Decimal numbers, +, -, *, /, and ordinary
-- Lua-style function calls all pass through completely unchanged --
-- nothing here re-implements arithmetic, only the small set of
-- conveniences (bracket styles, function names, log_N, **, e) that Lua's
-- own syntax doesn't already provide.
--
-- Euler's number `e` is supported two ways, both meaning e^(...): as a
-- function call, e(x+1), rewritten to Calc_one_var_exp(x+1); and as a bare
-- constant combined with **, e**{x+1} (normalized to e^(x+1) by the
-- bracket and ** steps), rewritten to a plain Calc_one_var_e^(x+1). The
-- bare-constant form is scanned by hand rather than gsub, specifically to
-- avoid corrupting Lua's own scientific-notation numeric literals
-- (5e3, 1.5e-10), which use the same letter for something unrelated.
--
-- "j" is the imaginary unit (sqrt(-1), engineering convention): a bare
-- constant like "e", but safe to gsub directly since "j" never appears in
-- Lua's own numeric-literal grammar. Complex arithmetic (+ - * / ^) works
-- once "j" is in play, and so do sin/cos/tan/asin/acos/atan/sinh/cosh/
-- tanh/exp/ln/sqrt/abs -- each dispatches to a Complex-aware formula
-- whenever its argument actually is Complex, and otherwise behaves
-- exactly as before (still returns a plain number). ln, sqrt, asin, and
-- acos additionally auto-promote an out-of-real-domain argument to
-- Complex (sqrt(-4) = 2i, ln(-1) = i*pi, asin(2) is complex) instead of
-- silently returning nan.
--
-- "root_N(...)" is the N-th (possibly decimal) root of (...), e.g.
-- root_3(x) or root_6.13(x) -- the same prefix-rewrite trick as log_N,
-- built on the new Root_any(index, radicand) function.
--
-- Whitespace anywhere in `expr` is insignificant -- stripped up front, so
-- "x**7+cos(x)" and "x**7   +   cos( x )" compile identically.
--
-- NOT yet implemented (deliberately, as a separate next step): the
-- inv(func(...)) wrapper for inverse functions (e.g. inv(tan(x)) as
-- shorthand for atan(x)) -- unlike log_N, that needs the matching
-- close-paren of func(...) to be found and removed, not just a prefix
-- rewrite. The bracket scanner above is most of what that needs; it's a
-- real next increment, not folded into this one.
--
-- Claude: returns TWO values -- f, and the normalized Lua source string
-- f was actually compiled from (e.g. "Calc_one_var_e^x+7" for "e**x + 7",
-- vs "Calc_one_var_e^(x+7)" for "e**{x + 7}"). Operator precedence here is
-- fully deterministic (exponentiation binds tighter than +/-, exactly like
-- standard math notation), so there is no real parsing ambiguity to warn
-- about -- but a user translating from handwritten superscripts, where
-- there's no visual boundary on where an exponent ends, can still type
-- something other than what they meant. The second return value is the
-- honest fix for that: it lets a caller SEE exactly how grouping resolved
-- (present or absent parentheses) before trusting any numeric output,
-- rather than a heuristic warning that would have to fire on every
-- ordinary polynomial (x**2 + 1 is not a mistake). A caller that only
-- wants f can just ignore the second value, same as always.
local function Compile_Expression(expr)
    Core.assert_nonblank_string(expr, "Compile_Expression: expr")
    local normalized = strip_whitespace(expr)
    normalized = normalize_brackets(normalized)
    normalized = translate_log_base(normalized)
    normalized = translate_root_index(normalized)
    normalized = translate_function_names(normalized)
    normalized = translate_e_call(normalized)
    normalized = translate_bare_e(normalized)
    normalized = translate_imaginary_unit(normalized)
    normalized = normalized:gsub("%*%*", "^")
    assert_allowed_identifiers(normalized, expr)

    local source = "return function(x) return " .. normalized .. " end"
    local chunk, load_err = load(source, "Compile_Expression", "t", COMPILE_ENV)
    assert(chunk, "Compile_Expression: could not compile \"" .. expr .. "\" -- " ..
        tostring(load_err))
    local ok, f = pcall(chunk)
    assert(ok, "Compile_Expression: error while building the function from \"" ..
        expr .. "\" -- " .. tostring(f))
    return f, normalized
end

-- Claude: accepts EITHER an already-compiled Lua function (e.g. the
-- first return value of Compile_expression itself, or any function a
-- caller built some other way) OR a symbolic-expression string like
-- "sin(x**2)", compiling the latter automatically. This is what lets
-- dif_numeric.lua/integration_numeric.lua's public functions
-- (Derivative_at_point, Integral_definite, ...) accept the friendly
-- syntax directly, e.g. Derivative_at_point("sin(x**2)", 0.5), without
-- forcing every caller to remember to call Compile_expression first --
-- while still accepting a plain function unchanged, for a caller who
-- already has one.
local function Resolve_function(f, label)
    if type(f) == "string" then
        return (Compile_Expression(f))
    end
    Core.assert_callable(f, label)
    return f
end

return {
    Compile_expression = Compile_Expression,
    Root_any = Root_any,
    Resolve_function = Resolve_function,
}
