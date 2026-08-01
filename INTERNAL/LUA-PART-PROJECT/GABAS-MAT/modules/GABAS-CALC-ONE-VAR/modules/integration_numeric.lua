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
-- section 50). Infinite intervals, endpoint singularities (tanh-sinh),
-- Cauchy principal values, oscillatory-specific methods, and complex
-- integration contours are all explicitly NOT implemented here -- each is
-- its own specialized method in the specification (sections 12-14), not a
-- mode of the ordinary adaptive rule, and folding them in now would risk
-- exactly what the document warns against: silently returning a wrong
-- answer for a case this rule was never meant to handle. Same status as
-- inv(): acknowledged, deferred, not guessed at.
local Complex = require("GABAS-LINAL.modules.complex")
local Core = require("GABAS-CALC-ONE-VAR.modules.core")
local Quadrature = require("GABAS-CALC-ONE-VAR.modules.quadrature")

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
    Core.assert_callable(f, "Integral_definite: f")
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

return {
    Integral_definite = Integral_definite,
}
