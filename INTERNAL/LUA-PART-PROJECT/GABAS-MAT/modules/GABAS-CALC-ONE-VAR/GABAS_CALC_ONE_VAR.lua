-- [Calc One Var]
-- One-variable calculus. Starts with Compile_Expression: the prerequisite
-- for everything else here (numeric derivatives today, symbolic ones
-- later) -- it turns a string the user typed in a friendly, calculator-
-- like syntax into a real, callable Lua function of one variable x, since
-- every numerical method (finite differences, root finding, quadrature...)
-- needs an actual function to call, not a piece of text.
local Complex = require("GABAS-LINAL.modules.complex")

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
-- as an ordinary function call.
local function translate_e_call(expr)
    return (expr:gsub("%f[%a]e%(", "math.exp("))
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
local FUNCTION_ALIASES = {
    sin = "math.sin", cos = "math.cos", tan = "math.tan",
    asin = "math.asin", acos = "math.acos", atan = "math.atan",
    sinh = "Calc_one_var_sinh", cosh = "Calc_one_var_cosh", tanh = "Calc_one_var_tanh",
    abs = "math.abs", sqrt = "math.sqrt", exp = "math.exp", ln = "math.log",
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
-- Complex arithmetic (+, -, *, /, ^) works through Complex's operator
-- metamethods once this constant is in play; math.* functions (sin, exp,
-- ...) do NOT accept a Complex argument -- that is a separate, not-yet-
-- built extension, same status as inv().
local function translate_imaginary_unit(expr)
    return (expr:gsub("%f[%a]j%f[%A]", "Calc_one_var_j"))
end

-- Claude: sinh/cosh/tanh from their exponential definitions, NOT Lua's
-- own math.sinh/cosh/tanh -- those were made optional (LUA_COMPAT_5_2) as
-- of Lua 5.3 and aren't guaranteed to exist in every build, including
-- whatever Lua LuaTeX embeds, which is the actual target for this module.
-- Defining these ourselves means this never depends on a compile flag we
-- don't control.
local function Calc_one_var_sinh(x) return (math.exp(x) - math.exp(-x)) / 2 end
local function Calc_one_var_cosh(x) return (math.exp(x) + math.exp(-x)) / 2 end
local function Calc_one_var_tanh(x) return Calc_one_var_sinh(x) / Calc_one_var_cosh(x) end
local function Calc_one_var_log_base(base, x) return math.log(x, base) end

-- The N-th (possibly non-integer) root of radicand, real-valued whenever
-- one exists: directly as radicand^(1/index) when radicand >= 0, or, for a
-- negative radicand, only when index is an odd integer (the one case a
-- real root of a negative number is still well-defined, e.g. the cube
-- root of -8 is -2) -- computed as -((-radicand)^(1/index)) to avoid
-- raising a negative base to a fractional power, which Lua's `^` cannot
-- do (it would silently produce nan).
local function Root_any(index, radicand)
    assert(type(index) == "number" and index == index, "Root_any: index must be a real number.")
    assert(index ~= 0, "Root_any: index cannot be zero.")
    assert(type(radicand) == "number" and radicand == radicand, "Root_any: radicand must be a real number.")
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
-- NOT Lua's real global table (_G). Compiled user expressions only ever
-- see `math` and the handful of helpers above; they can't read or write
-- anything else in the process, and nothing in the process can be
-- shadowed by whatever a user happens to type.
local COMPILE_ENV = {
    math = math,
    Calc_one_var_log_base = Calc_one_var_log_base,
    Calc_one_var_sinh = Calc_one_var_sinh,
    Calc_one_var_cosh = Calc_one_var_cosh,
    Calc_one_var_tanh = Calc_one_var_tanh,
    Calc_one_var_e = Calc_one_var_e,
    Calc_one_var_j = Calc_one_var_j,
    Root_any = Root_any,
}

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
-- function call, e(x+1), rewritten to math.exp(x+1); and as a bare
-- constant combined with **, e**{x+1} (normalized to e^(x+1) by the
-- bracket and ** steps), rewritten to a plain Calc_one_var_e^(x+1). The
-- bare-constant form is scanned by hand rather than gsub, specifically to
-- avoid corrupting Lua's own scientific-notation numeric literals
-- (5e3, 1.5e-10), which use the same letter for something unrelated.
--
-- "j" is the imaginary unit (sqrt(-1), engineering convention): a bare
-- constant like "e", but safe to gsub directly since "j" never appears in
-- Lua's own numeric-literal grammar. Complex arithmetic (+ - * / ^) works
-- once "j" is in play; math.* functions do not accept a Complex argument.
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
    assert(type(expr) == "string" and expr:match("%S"),
        "Compile_Expression: expr must be a non-blank string.")
    local normalized = strip_whitespace(expr)
    normalized = normalize_brackets(normalized)
    normalized = translate_log_base(normalized)
    normalized = translate_root_index(normalized)
    normalized = translate_function_names(normalized)
    normalized = translate_e_call(normalized)
    normalized = translate_bare_e(normalized)
    normalized = translate_imaginary_unit(normalized)
    normalized = normalized:gsub("%*%*", "^")

    local source = "return function(x) return " .. normalized .. " end"
    local chunk, load_err = load(source, "Compile_Expression", "t", COMPILE_ENV)
    assert(chunk, "Compile_Expression: could not compile \"" .. expr .. "\" -- " ..
        tostring(load_err))
    local ok, f = pcall(chunk)
    assert(ok, "Compile_Expression: error while building the function from \"" ..
        expr .. "\" -- " .. tostring(f))
    return f, normalized
end

return {
    Compile_expression = Compile_Expression,
    Root_any = Root_any,
}
