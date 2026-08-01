-- [Integration Numeric]
-- Integral_definite(f, a, b, opts) -> value, estimated_error, status,
-- diagnostics
--
-- Adaptive Gauss-Kronrod 21/10 quadrature (modules/quadrature.lua),
-- following the technical document's own recommended architecture:
-- repeatedly bisect the subinterval with the largest estimated error
-- until the global (compensated) error estimate meets
-- max(abs_tol, rel_tol * |integral|).
--
-- Claude: scoped deliberately to REAL, FINITE limits and a real-or-
-- Complex-valued integrand -- exactly the "ordinary case" the document
-- treats as the primary engine (its own final recommended architecture,
-- section 50). Endpoint singularities (tanh-sinh), Cauchy principal
-- values, oscillatory-specific methods, and complex integration contours
-- are all explicitly NOT implemented here -- each is its own specialized
-- method in the specification (sections 12-14), not a mode of the
-- ordinary adaptive rule, and folding them in now would risk exactly what
-- the document warns against: silently returning a wrong answer for a
-- case this rule was never meant to handle. Same status as inv():
-- acknowledged, deferred, not guessed at.
--
-- Infinite intervals ARE handled, by Integral_infinite below -- but its
-- own algebraic transform has a known, tested limitation worth flagging
-- here too: a SLOWLY (power-law) decaying tail turns into a genuine
-- endpoint singularity in the transformed variable (verified: integral
-- of x^-p from 1 to infinity transforms to integral of (1-t)^(p-2) on
-- [0,1], singular at t=1 whenever p<2), which this ordinary adaptive
-- rule -- with no singularity-specific handling -- cannot always resolve
-- as p approaches the convergence boundary (p=1). It correctly reports a
-- non-SUCCESS status rather than a false convergence in that case (tested
-- down to p=1.5); a genuinely fast-decaying tail (exponential, Gaussian,
-- 1/(1+x^2), ...) converges cleanly. tanh-sinh, once built, is the
-- documented fallback for exactly this combination.
local Complex = require("GABAS-LINAL.modules.complex")
local Core = require("GABAS-CALC-ONE-VAR.modules.core")
local Quadrature = require("GABAS-CALC-ONE-VAR.modules.quadrature")
local TanhSinh = require("GABAS-CALC-ONE-VAR.modules.tanh_sinh")
local CompileExpression = require("GABAS-CALC-ONE-VAR.modules.compile_expression")

-- Claude: compensated (Neumaier/Kahan) summation over a list of plain
-- reals or Complex values, real and imaginary parts compensated
-- separately (technical document, section 10) -- used for the GLOBAL
-- running integral/error total, recomputed fresh from every active
-- interval on each convergence check rather than updated incrementally
-- (subtract-old-add-new), avoiding exactly the accumulated rounding the
-- document warns about: "naive summation... can lose small contributions
-- when the current sum is much larger."
local function kahan_sum(values)
    local sum_re, comp_re = 0, 0
    local sum_im, comp_im = 0, 0
    local any_complex = false
    for i = 1, #values do
        local v = values[i]
        local re, im
        if Complex.Is_complex(v) then
            re, im = v.re, v.im
            any_complex = true
        else
            re, im = v, 0
        end
        local y = re - comp_re
        local t = sum_re + y
        comp_re = (t - sum_re) - y
        sum_re = t

        local yi = im - comp_im
        local ti = sum_im + yi
        comp_im = (ti - sum_im) - yi
        sum_im = ti
    end
    if any_complex then
        return Complex.Complex(sum_re, sum_im)
    end
    return sum_re
end

local function kahan_sum_real(values)
    local sum, comp = 0, 0
    for i = 1, #values do
        local y = values[i] - comp
        local t = sum + y
        comp = (t - sum) - y
        sum = t
    end
    return sum
end

-- opts (all optional): abs_tol (default 1e-10), rel_tol (default 1e-8),
-- max_intervals (default 100).
--
-- Finding the worst-error interval each iteration is a plain linear scan,
-- not a heap -- deliberately: one interval refinement already costs 21+
-- integrand evaluations, which dominates the total cost long before a
-- linear scan over (at most max_intervals) entries would matter, and the
-- simpler structure is less likely to hide a bookkeeping bug than a
-- hand-rolled heap in a first implementation.
local function Integral_definite(f, a, b, opts)
    opts = opts or {}
    assert(type(opts) == "table", "Integral_definite: opts must be a table.")
    f = CompileExpression.Resolve_function(f, "Integral_definite: f")
    Core.assert_finite_number(a, "Integral_definite: a")
    Core.assert_finite_number(b, "Integral_definite: b")

    local abs_tol = opts.abs_tol or 1e-10
    local rel_tol = opts.rel_tol or 1e-8
    Core.assert_nonneg_number(abs_tol, "Integral_definite: opts.abs_tol")
    Core.assert_nonneg_number(rel_tol, "Integral_definite: opts.rel_tol")
    assert(abs_tol > 0 or rel_tol > 0,
        "Integral_definite: at least one of opts.abs_tol/opts.rel_tol must be positive.")
    local max_intervals = opts.max_intervals or 100
    Core.assert_positive_integer(max_intervals, "Integral_definite: opts.max_intervals")

    -- Zero-length interval: the integral is 0 by definition, without ever
    -- evaluating f -- important when f is undefined exactly at a=b.
    if a == b then
        return 0, 0, "SUCCESS", {evaluations = 0, intervals = 0}
    end

    local orientation = 1
    if b < a then
        a, b = b, a
        orientation = -1
    end

    local intervals = {}
    local evaluations = 0

    local function refine(left, right)
        local result, error, resabs, resasc, evals, err_msg =
            Quadrature.Local_Gauss_Kronrod_21(f, left, right)
        evaluations = evaluations + evals
        if result == nil then
            return nil, err_msg
        end
        return {left = left, right = right, integral = result, error = error}
    end

    local first, first_err = refine(a, b)
    if first == nil then
        return nil, nil, "INTEGRAND_EVALUATION_FAILURE",
            {evaluations = evaluations, intervals = 0, message = first_err}
    end
    intervals[1] = first

    local function global_totals()
        local sub_integrals, sub_errors = {}, {}
        for i = 1, #intervals do
            sub_integrals[i] = intervals[i].integral
            sub_errors[i] = intervals[i].error
        end
        return kahan_sum(sub_integrals), kahan_sum_real(sub_errors)
    end

    local global_integral, global_error = global_totals()
    local status = "SUCCESS"

    while true do
        local target = math.max(abs_tol, rel_tol * Quadrature.magnitude(global_integral))
        if global_error <= target then
            break
        end
        if #intervals >= max_intervals then
            status = "INTERVAL_LIMIT_REACHED"
            break
        end

        local worst_idx, worst_err = 1, intervals[1].error
        for i = 2, #intervals do
            if intervals[i].error > worst_err then
                worst_idx, worst_err = i, intervals[i].error
            end
        end
        local worst = intervals[worst_idx]
        local mid = worst.left + (worst.right - worst.left) / 2
        if mid == worst.left or mid == worst.right then
            status = "FLOATING_POINT_INTERVAL_COLLAPSE"
            break
        end

        local left_child = refine(worst.left, mid)
        if left_child == nil then
            status = "INTEGRAND_EVALUATION_FAILURE"
            break
        end
        local right_child = refine(mid, worst.right)
        if right_child == nil then
            status = "INTEGRAND_EVALUATION_FAILURE"
            break
        end

        intervals[worst_idx] = left_child
        intervals[#intervals + 1] = right_child
        global_integral, global_error = global_totals()
    end

    return orientation * global_integral, global_error, status,
        {evaluations = evaluations, intervals = #intervals}
end

-- Integral_infinite(f, a, b, opts) -> value, estimated_error, status,
-- diagnostics
--
-- Definite integral where either or both limits may be infinite
-- (math.huge / -math.huge), via the algebraic transformations of section
-- 13: a semi-infinite bound uses x = a + t/(1-t) (t in [0,1), dx/dt =
-- 1/(1-t)^2); the whole real line uses x = t/(1-t^2) (t in (-1,1), dx/dt
-- = (1+t^2)/(1-t^2)^2). Both transforms are then handed to the SAME
-- adaptive Integral_definite above -- Gauss-Kronrod nodes never land
-- exactly on an interval's endpoint (the largest Kronrod abscissa is
-- ~0.9957, not 1), so the transformed integrand is never evaluated
-- exactly at t=+-1, where the Jacobian would overflow or the map would
-- reach the literal infinite limit.
--
-- If both limits happen to already be finite, this delegates to
-- Integral_definite directly (no transform needed) -- safe to call
-- whenever a bound MIGHT be infinite, without checking first.
--
-- Claude: requires a < b (including for infinite bounds) -- this also
-- rejects the "two identical infinite endpoints" case the specification
-- explicitly calls out (section 18.5) as an input error, since
-- math.huge < math.huge is false.
local function Integral_infinite(f, a, b, opts)
    f = CompileExpression.Resolve_function(f, "Integral_infinite: f")
    assert(type(a) == "number" and a == a, "Integral_infinite: a must be a number (possibly +-math.huge).")
    assert(type(b) == "number" and b == b, "Integral_infinite: b must be a number (possibly +-math.huge).")
    assert(a < b, "Integral_infinite: a must be less than b.")

    local a_finite = a ~= -math.huge and a ~= math.huge
    local b_finite = b ~= -math.huge and b ~= math.huge

    if a_finite and b_finite then
        return Integral_definite(f, a, b, opts)
    end

    if a_finite and b == math.huge then
        local g = function(t)
            local one_minus_t = 1 - t
            return f(a + t / one_minus_t) / (one_minus_t * one_minus_t)
        end
        return Integral_definite(g, 0, 1, opts)
    end

    if a == -math.huge and b_finite then
        local g = function(t)
            local one_minus_t = 1 - t
            return f(b - t / one_minus_t) / (one_minus_t * one_minus_t)
        end
        return Integral_definite(g, 0, 1, opts)
    end

    local g = function(t)
        local one_minus_t2 = 1 - t * t
        return f(t / one_minus_t2) * (1 + t * t) / (one_minus_t2 * one_minus_t2)
    end
    return Integral_definite(g, -1, 1, opts)
end

-- Integral_principal_value(f, a, b, c, opts) -> value, estimated_error,
-- status, diagnostics
--
-- The Cauchy principal value of the integral of f from a to b, where f
-- has a singularity exactly at the interior point c (a < c < b):
-- PV = lim(eps->0) [integral from a to c-eps + integral from c+eps to b]
-- (technical document, section 2.4).
--
-- Claude: the technique here is symmetric folding, substituting t = x-c
-- and defining phi(t) = f(c+t) + f(c-t) for t > 0. For a genuine simple-
-- pole singularity (f(x) ~ k/(x-c) near c, the classical PV case), the
-- divergent 1/t parts of f(c+t) and f(c-t) are exactly opposite and
-- cancel in the sum -- phi is then an ordinary, non-singular function of
-- t on (0, h], h = min(c-a, b-c), safe to hand straight to the existing
-- adaptive Integral_definite (whose Gauss-Kronrod nodes never land
-- exactly on t=0 anyway, so phi is never evaluated at the fold point
-- itself). Any leftover, non-symmetric portion of [a,b] outside
-- [c-h, c+h] is a completely ordinary integral with no singularity in it
-- at all, added on separately.
--
-- Verified against PV(1/x) on [-1,1] = 0, PV(1/x) on [-1,2] = ln(2)
-- (symmetric part folds to exactly 0, leaving only the ordinary [1,2]
-- tail), and PV((x+1)/x) on [-1,1] = 2. This is specifically a technique
-- for an ODD-type (simple-pole-like) singularity -- it is not a general
-- antidote for an arbitrary discontinuity, and does not (cannot) fold
-- away an EVEN-type singularity like 1/x^2, which genuinely has no
-- principal value; Integral_definite on the symmetric part correctly
-- reports failure/non-convergence for that case rather than a false one.
local function Integral_principal_value(f, a, b, c, opts)
    f = CompileExpression.Resolve_function(f, "Integral_principal_value: f")
    Core.assert_finite_number(a, "Integral_principal_value: a")
    Core.assert_finite_number(b, "Integral_principal_value: b")
    Core.assert_finite_number(c, "Integral_principal_value: c")
    assert(a < c and c < b, "Integral_principal_value: c must lie strictly between a and b.")

    local h = math.min(c - a, b - c)
    local phi = function(t) return f(c + t) + f(c - t) end

    local value, error, status, diag = Integral_definite(phi, 0, h, opts)
    local evaluations, intervals = diag.evaluations, diag.intervals

    if c - h > a then
        local lv, le, ls, ld = Integral_definite(f, a, c - h, opts)
        value = value + lv
        error = error + le
        evaluations = evaluations + ld.evaluations
        intervals = intervals + ld.intervals
        if ls ~= "SUCCESS" then status = ls end
    end
    if c + h < b then
        local rv, re, rs, rd = Integral_definite(f, c + h, b, opts)
        value = value + rv
        error = error + re
        evaluations = evaluations + rd.evaluations
        intervals = intervals + rd.intervals
        if rs ~= "SUCCESS" then status = rs end
    end

    return value, error, status, {evaluations = evaluations, intervals = intervals}
end

-- Claude: one trapezoidal sum at step h, over the DE-transformed
-- variable t -- see tanh_sinh.lua's De_map header for exactly which
-- nodes are skipped (and why that is always safe). lenient_eval is
-- responsible for turning a NaN/infinite/failed evaluation into a
-- harmless 0 contribution -- by construction this should essentially
-- never trigger for a genuinely integrable singularity (the
-- doubly-exponential weight already suppresses the integrand to
-- negligible before floating-point precision could overflow it), so
-- repeated triggering here is itself a signal the singularity is not
-- integrable, which the level-doubling convergence check in
-- Integral_tanh_sinh will correctly fail to satisfy.
local function tanh_sinh_level(lenient_eval, a, b, mid, half_range, h)
    local sum = select(2, TanhSinh.De_map(0)) * lenient_eval(mid)
    local evaluations = 1
    local k = 1
    local MAX_T = 12
    while true do
        local t = k * h
        if t > MAX_T then break end
        local delta, w = TanhSinh.De_map(t)
        if not delta then break end
        -- delta depends only on |t|, so this single De_map call serves
        -- both the node approaching b (t>0 side) and the node
        -- approaching a (t<0 side) -- see tanh_sinh.lua's header.
        local x_pos = b - half_range * delta
        local x_neg = a + half_range * delta
        sum = sum + w * (lenient_eval(x_pos) + lenient_eval(x_neg))
        evaluations = evaluations + 2
        k = k + 1
    end
    return sum * h * half_range, evaluations
end

-- Integral_tanh_sinh(f, a, b, opts) -> value, estimated_error, status,
-- diagnostics
--
-- Definite integral over a FINITE [a,b], via double-exponential (tanh-
-- sinh) quadrature (modules/tanh_sinh.lua) -- the specification's own
-- recommended fallback for an integrand with an INTEGRABLE algebraic or
-- logarithmic singularity at (or near) an endpoint, which the ordinary
-- adaptive Gauss-Kronrod rule (Integral_definite) cannot reliably
-- resolve (see Integral_infinite's own documented boundary case above,
-- which is exactly this problem in disguise). The doubly-exponential
-- decay of the DE weight overwhelms an algebraic blow-up in the
-- integrand near the endpoints, letting the trapezoidal rule -- applied
-- in the TRANSFORMED variable -- converge even though the ORIGINAL
-- integrand is unbounded there.
--
-- Claude: convergence is checked by level doubling (halving the step h
-- and comparing successive trapezoidal estimates), the classical
-- tanh-sinh termination strategy -- NOT the Gauss-Kronrod error
-- estimator, which assumes a different, polynomial-exactness-based
-- error model that does not apply to a trapezoidal rule on a
-- transformed variable. opts.tol (default 1e-10) is the convergence
-- threshold between consecutive levels; opts.max_levels (default 12,
-- i.e. h down to 1/2^12) bounds the work.
--
-- The interval's own center is checked strictly (any failure there
-- aborts immediately with a clear error) -- it is never expected to be
-- near a singularity, so a failure there signals something is
-- fundamentally wrong with f, not merely an endpoint difficulty.
--
-- Claude: this handles an ENDPOINT singularity (at a or at b) --
-- verified exactly (to machine precision) against integral of x^-0.5 on
-- [0,1] = 2 and integral of ln(x) on [0,1] = -1. An INTERIOR singularity
-- (e.g. |x|^-0.5 on [-1,1], singular at the interior point x=0) is a
-- DIFFERENT problem this function does not detect or handle on its own
-- -- split the interval at the singular point first and call this twice
-- (verified: splitting |x|^-0.5 at x=0 and summing both halves gives
-- exactly 4, matching the closed form). As the singularity gets close
-- enough to non-integrable (e.g. x^alpha near the alpha=-1 boundary),
-- there is a genuine, expected precision limit -- more levels stop
-- helping once the trapezoidal sum's own truncation error, not the
-- tolerance, is the binding constraint (verified: alpha=-0.99 plateaus
-- around a ~0.08% error regardless of opts.max_levels, correctly
-- reporting MAX_LEVELS_REACHED rather than a false SUCCESS).
local function Integral_tanh_sinh(f, a, b, opts)
    opts = opts or {}
    assert(type(opts) == "table", "Integral_tanh_sinh: opts must be a table.")
    f = CompileExpression.Resolve_function(f, "Integral_tanh_sinh: f")
    Core.assert_finite_number(a, "Integral_tanh_sinh: a")
    Core.assert_finite_number(b, "Integral_tanh_sinh: b")
    assert(a < b, "Integral_tanh_sinh: a must be less than b.")
    local tol = opts.tol or 1e-10
    Core.assert_positive_number(tol, "Integral_tanh_sinh: opts.tol")
    local max_levels = opts.max_levels or 12
    Core.assert_positive_integer(max_levels, "Integral_tanh_sinh: opts.max_levels")

    local mid = (a + b) / 2
    local half_range = (b - a) / 2

    local function lenient_eval(x)
        local ok, v = pcall(f, x)
        if not ok then return 0 end
        if Complex.Is_complex(v) then
            if v.re ~= v.re or v.im ~= v.im or v.re == math.huge or v.re == -math.huge
                or v.im == math.huge or v.im == -math.huge then
                return 0
            end
            return v
        end
        if type(v) ~= "number" or v ~= v or v == math.huge or v == -math.huge then
            return 0
        end
        return v
    end

    local center_ok, center_err = pcall(f, mid)
    assert(center_ok, "Integral_tanh_sinh: f(" .. tostring(mid) ..
        ") (the interval's center) failed to evaluate -- " .. tostring(center_err))

    local total_evaluations = 0
    local h = 1
    local estimate = nil
    for level = 1, max_levels do
        local level_estimate, evals = tanh_sinh_level(lenient_eval, a, b, mid, half_range, h)
        total_evaluations = total_evaluations + evals
        if estimate ~= nil then
            local diff = Quadrature.magnitude(level_estimate - estimate)
            if diff <= tol then
                return level_estimate, diff, "SUCCESS", {evaluations = total_evaluations, levels = level}
            end
        end
        estimate = level_estimate
        h = h / 2
    end

    return estimate, nil, "MAX_LEVELS_REACHED", {evaluations = total_evaluations, levels = max_levels}
end

return {
    Integral_definite = Integral_definite,
    Integral_infinite = Integral_infinite,
    Integral_principal_value = Integral_principal_value,
    Integral_tanh_sinh = Integral_tanh_sinh,
}
