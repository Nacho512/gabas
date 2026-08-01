-- [Optimization Numeric]
-- Numerical optimization for a REAL one-variable function: finding an x
-- where f(x) is smallest or largest, NOT where f(x) = 0 (that is
-- roots_numeric's job). Two families:
--
--   Golden_section_minimize/maximize -- derivative-free bracketing,
--   assumes f is unimodal on [a,b] (exactly one interior extremum).
--
--   Local_minimum/Local_maximum -- derivative-based: an extremum of a
--   differentiable f occurs where f'(x) = 0, so finding one is root-
--   finding on f' -- this reuses Derivative_at_point (exact, via Dual
--   numbers) and roots_numeric.Secant directly, with no new derivative
--   machinery of its own. Classifies what it finds by checking the sign
--   of f' just to either side (min: -to-+, max: +to--) rather than
--   assuming a critical point IS the kind of extremum the caller asked
--   for -- an inflection point (f'=0 but no sign change, e.g. x^3 at 0)
--   is reported honestly as "not a minimum"/"not a maximum", never
--   silently mislabeled.
--
--   Global_minimum/Global_maximum -- on a CLOSED interval, without
--   assuming unimodality: samples f' across the interval, refines every
--   detected sign change into a precise critical point (reusing the SAME
--   Secant-on-f' building block Local_minimum uses), then compares f at
--   every critical point found AND the two endpoints (the extreme value
--   theorem, in code). Honest limitation, stated plainly rather than
--   hidden: a coarse enough sample grid can still miss a very narrow
--   feature entirely between two samples -- opts.samples raises the
--   resolution at the cost of more evaluations, but this is a sampling
--   method, not a certified one (that would need interval arithmetic, a
--   different, future technique).
local Core = require("GABAS-CALC-ONE-VAR.modules.core")
local CompileExpression = require("GABAS-CALC-ONE-VAR.modules.compile_expression")
local DifNumeric = require("GABAS-CALC-ONE-VAR.modules.dif_numeric")
local RootsNumeric = require("GABAS-CALC-ONE-VAR.modules.roots_numeric")

local GOLDEN_RATIO_CONJUGATE = (math.sqrt(5) - 1) / 2

-- Golden_section_minimize(f, a, b, opts) -> x_min, f_min, status,
-- diagnostics
--
-- Classical golden-section search: narrows the bracket [a,b] by
-- comparing f at two interior points placed a golden-ratio fraction in
-- from each end, discarding whichever side cannot contain the minimum.
-- Only ONE new evaluation is needed per iteration -- the other interior
-- point is always reused from the previous step, the well-known
-- self-similarity property of this exact ratio. No derivative is used
-- at all, so this works for a non-differentiable (but still unimodal)
-- f, unlike Local_minimum below.
--
-- ASSUMES f is unimodal on [a,b] (has exactly one interior minimum,
-- decreasing then increasing) -- not verified (verifying unimodality in
-- general is itself a hard problem), just assumed, same as every other
-- specialized method in this package that requires a specific structural
-- property of its input.
--
-- opts (all optional): tol (default 1e-10, the final bracket width),
-- max_iterations (default 200 -- golden-section shrinks the bracket by
-- a factor of ~0.618 per iteration, so this comfortably allows shrinking
-- to far below any reasonable tolerance).
--
-- Claude: the bracket-width stopping criterion itself is exact, but the
-- achievable PRECISION IN x is limited by how much f actually varies
-- near the optimum relative to f's own absolute scale there -- a large
-- additive constant in f (e.g. "-(x-2)^2 + 1000") makes the quadratic
-- term near the optimum drop below that constant's own floating-point
-- ULP, so f(c) and f(d) can round to the identical value for a whole
-- neighborhood of x, starving the comparison of any real signal.
-- Verified directly (see the test suite): this is a genuine
-- floating-point precision limit of any comparison-based optimizer, not
-- specific to this implementation.
local function Golden_section_minimize(f, a, b, opts)
    opts = opts or {}
    assert(type(opts) == "table", "Golden_section_minimize: opts must be a table.")
    f = CompileExpression.Resolve_function(f, "Golden_section_minimize: f")
    Core.assert_finite_number(a, "Golden_section_minimize: a")
    Core.assert_finite_number(b, "Golden_section_minimize: b")
    assert(a ~= b, "Golden_section_minimize: a and b must be distinct.")
    if b < a then a, b = b, a end

    local tol = opts.tol or 1e-10
    Core.assert_positive_number(tol, "Golden_section_minimize: opts.tol")
    local max_iterations = opts.max_iterations or 200
    Core.assert_positive_integer(max_iterations, "Golden_section_minimize: opts.max_iterations")

    local c = b - GOLDEN_RATIO_CONJUGATE * (b - a)
    local d = a + GOLDEN_RATIO_CONJUGATE * (b - a)
    local fc, fd = f(c), f(d)
    Core.assert_finite_number(fc, "Golden_section_minimize: f(c)")
    Core.assert_finite_number(fd, "Golden_section_minimize: f(d)")

    local status = "SUCCESS"
    local iteration = 0
    while (b - a) > tol do
        iteration = iteration + 1
        if iteration > max_iterations then
            status = "MAX_ITERATIONS_REACHED"
            break
        end
        if fc < fd then
            b, d, fd = d, c, fc
            c = b - GOLDEN_RATIO_CONJUGATE * (b - a)
            fc = f(c)
            Core.assert_finite_number(fc, "Golden_section_minimize: f(c)")
        else
            a, c, fc = c, d, fd
            d = a + GOLDEN_RATIO_CONJUGATE * (b - a)
            fd = f(d)
            Core.assert_finite_number(fd, "Golden_section_minimize: f(d)")
        end
    end

    local x_min, f_min = c, fc
    if fd < f_min then x_min, f_min = d, fd end
    return x_min, f_min, status, {iterations = iteration, bracket_width = b - a}
end

-- Golden_section_maximize(f, a, b, opts) -> x_max, f_max, status,
-- diagnostics
--
-- Maximizing f is minimizing -f -- delegates to Golden_section_minimize
-- on the negated function rather than duplicating the search logic, and
-- negates the returned value back.
local function Golden_section_maximize(f, a, b, opts)
    f = CompileExpression.Resolve_function(f, "Golden_section_maximize: f")
    local negated = function(x) return -f(x) end
    local x_max, neg_f_max, status, diag = Golden_section_minimize(negated, a, b, opts)
    return x_max, -neg_f_max, status, diag
end

-- Claude: g(x) = f'(x), computed exactly via Derivative_at_point (Dual
-- numbers) rather than a finite difference -- the shared building block
-- both find_critical_point_and_classify and the global-search sampler
-- below reduce to. f must therefore be Dual-compatible (a symbolic-
-- expression string, auto-compiled, or an already-compiled function
-- built only from Compile_expression's own Dual-aware primitives) --
-- same requirement Newton in roots_numeric.lua already has, for the
-- same reason.
local function derivative_fn(f)
    return function(x) return (DifNumeric.Derivative_at_point(f, x)) end
end

-- Claude: finds ONE critical point near x0 (root-finding on f' via
-- Secant, seeded with x0 and a nearby second point) and classifies it by
-- checking the sign of f' immediately to either side -- MINIMUM
-- (- to +), MAXIMUM (+ to -), or NEITHER (an inflection point like x^3
-- at 0, where f'=0 but does not change sign). The probe distance for
-- that sign check uses the same cube-root-of-machine-epsilon scaling as
-- Verify_derivative_at_point in dif_numeric.lua -- small enough to probe
-- genuinely local behavior, large enough to stay clear of floating-point
-- noise.
local function find_critical_point_and_classify(f, x0, opts)
    local g = derivative_fn(f)
    local step = opts.initial_step or 0.01 * math.max(1, math.abs(x0))
    local x_star, _, status, diag = RootsNumeric.Secant(g, x0, x0 + step, opts)
    if status ~= "SUCCESS" then
        return x_star, nil, nil, status, diag
    end
    local eps = (Core.MACHINE_EPSILON ^ (1 / 3)) * math.max(1, math.abs(x_star))
    local slope_left = g(x_star - eps)
    local slope_right = g(x_star + eps)
    local kind = "NEITHER"
    if slope_left < 0 and slope_right > 0 then
        kind = "MINIMUM"
    elseif slope_left > 0 and slope_right < 0 then
        kind = "MAXIMUM"
    end
    return x_star, f(x_star), kind, "SUCCESS", diag
end

-- Local_minimum(f, x0, opts) -> x_min, f_min, status, diagnostics
--
-- The local minimum of f nearest to the starting guess x0, found by
-- root-finding on f' (see find_critical_point_and_classify above) and
-- verifying the result genuinely IS a minimum, not just any critical
-- point. opts are passed straight through to the underlying Secant call
-- (opts.tol, opts.max_iterations), plus opts.initial_step (default
-- 0.01*max(1,|x0|)) for Secant's own second starting point.
--
-- status is "SUCCESS" only when a genuine local minimum was found;
-- "NOT_A_MINIMUM" if the critical point found is actually a maximum or
-- an inflection point (never silently mislabeled); otherwise whatever
-- honest failure status the underlying Secant search reported
-- (STALLED, DIVERGED, MAX_ITERATIONS_REACHED).
local function Local_minimum(f, x0, opts)
    opts = opts or {}
    assert(type(opts) == "table", "Local_minimum: opts must be a table.")
    f = CompileExpression.Resolve_function(f, "Local_minimum: f")
    Core.assert_finite_number(x0, "Local_minimum: x0")
    local x_star, f_star, kind, status, diag = find_critical_point_and_classify(f, x0, opts)
    if status ~= "SUCCESS" then
        return x_star, nil, status, diag
    end
    if kind ~= "MINIMUM" then
        return x_star, f_star, "NOT_A_MINIMUM", diag
    end
    return x_star, f_star, "SUCCESS", diag
end

-- Local_maximum(f, x0, opts) -> x_max, f_max, status, diagnostics
--
-- The mirror image of Local_minimum -- same critical-point search
-- (a maximum of f is also where f'=0), but requires the sign pattern to
-- classify as MAXIMUM rather than MINIMUM.
local function Local_maximum(f, x0, opts)
    opts = opts or {}
    assert(type(opts) == "table", "Local_maximum: opts must be a table.")
    f = CompileExpression.Resolve_function(f, "Local_maximum: f")
    Core.assert_finite_number(x0, "Local_maximum: x0")
    local x_star, f_star, kind, status, diag = find_critical_point_and_classify(f, x0, opts)
    if status ~= "SUCCESS" then
        return x_star, nil, status, diag
    end
    if kind ~= "MAXIMUM" then
        return x_star, f_star, "NOT_A_MAXIMUM", diag
    end
    return x_star, f_star, "SUCCESS", diag
end

-- Claude: samples f' at `samples` evenly spaced points across [a,b] and
-- refines every detected sign change (via the same Secant-on-f' call
-- find_critical_point_and_classify uses) into a precise critical x --
-- deliberately not reusing find_critical_point_and_classify itself,
-- since a global scan doesn't need its classification step (every
-- critical point found here gets compared by f-VALUE against every
-- other candidate below, not by kind). A sample that lands (to floating-
-- point precision) exactly on a critical point is also caught directly,
-- rather than relying on a sign change either side of it.
local function scan_critical_points(f, a, b, samples)
    local g = derivative_fn(f)
    local step = (b - a) / samples
    local candidates = {}
    local prev_x, prev_slope = a, g(a)
    if prev_slope == 0 then candidates[#candidates + 1] = a end
    for i = 1, samples do
        local x = a + i * step
        local slope = g(x)
        if slope == 0 then
            candidates[#candidates + 1] = x
        elseif (prev_slope < 0) ~= (slope < 0) then
            local root, _, status = RootsNumeric.Secant(g, prev_x, x)
            if status == "SUCCESS" then
                candidates[#candidates + 1] = root
            end
        end
        prev_x, prev_slope = x, slope
    end
    return candidates
end

-- Global_minimum(f, a, b, opts) -> x_min, f_min, status, diagnostics
--
-- The smallest value of f on the CLOSED interval [a,b] -- the extreme
-- value theorem, in code: every critical point found by scanning [a,b]
-- (see scan_critical_points above), plus the two endpoints themselves
-- (a continuous function's global extremum on a closed interval, if not
-- at a critical point, must be at an endpoint), compared by f-value.
-- Does NOT assume unimodality, unlike Golden_section_minimize.
--
-- opts (all optional): samples (default 50, the resolution of the
-- initial critical-point scan -- see this file's own header for the
-- honest limitation this implies).
local function Global_minimum(f, a, b, opts)
    opts = opts or {}
    assert(type(opts) == "table", "Global_minimum: opts must be a table.")
    f = CompileExpression.Resolve_function(f, "Global_minimum: f")
    Core.assert_finite_number(a, "Global_minimum: a")
    Core.assert_finite_number(b, "Global_minimum: b")
    assert(a < b, "Global_minimum: a must be less than b.")
    local samples = opts.samples or 50
    Core.assert_positive_integer(samples, "Global_minimum: opts.samples")

    local candidates = scan_critical_points(f, a, b, samples)
    candidates[#candidates + 1] = a
    candidates[#candidates + 1] = b

    local best_x, best_f = nil, nil
    for _, x in ipairs(candidates) do
        local fx = f(x)
        if best_f == nil or fx < best_f then
            best_x, best_f = x, fx
        end
    end
    return best_x, best_f, "SUCCESS", {candidates = #candidates, samples = samples}
end

-- Global_maximum(f, a, b, opts) -> x_max, f_max, status, diagnostics
--
-- Maximizing f is minimizing -f, same delegation as
-- Golden_section_maximize.
local function Global_maximum(f, a, b, opts)
    f = CompileExpression.Resolve_function(f, "Global_maximum: f")
    local negated = function(x) return -f(x) end
    local x_max, neg_f_max, status, diag = Global_minimum(negated, a, b, opts)
    return x_max, -neg_f_max, status, diag
end

return {
    Golden_section_minimize = Golden_section_minimize,
    Golden_section_maximize = Golden_section_maximize,
    Local_minimum = Local_minimum,
    Local_maximum = Local_maximum,
    Global_minimum = Global_minimum,
    Global_maximum = Global_maximum,
}
