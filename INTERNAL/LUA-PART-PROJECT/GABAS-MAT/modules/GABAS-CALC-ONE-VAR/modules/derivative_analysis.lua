-- [Derivative Analysis]
-- What the derivative reveals about f on a closed interval [a,b]: where
-- its critical points are and what kind each is, where f is increasing
-- or decreasing, and the tangent line / linear approximation at a point.
-- Built entirely on already-existing pieces -- modules/critical_points.lua
-- (shared with optimization_numeric.lua) for finding and classifying
-- critical points, and Derivative_at_point for the tangent line -- no
-- new derivative machinery of its own.
--
-- Claude: deliberately does NOT include concavity/inflection-point
-- analysis (concave-up/down intervals, inflection points) -- that needs
-- the SECOND derivative, which this package has no general way to
-- compute yet (Derivative_at_point's Dual numbers give an EXACT first
-- derivative; getting a second derivative the same way would need
-- nested Dual numbers, real new infrastructure, not something to bolt on
-- as a side effect of this file). Acknowledged, deferred, not guessed
-- at -- same status as the oscillatory-integral extension was before it
-- got built properly.
local Core = require("GABAS-CALC-ONE-VAR.modules.core")
local CompileExpression = require("GABAS-CALC-ONE-VAR.modules.compile_expression")
local DifNumeric = require("GABAS-CALC-ONE-VAR.modules.dif_numeric")
local CriticalPoints = require("GABAS-CALC-ONE-VAR.modules.critical_points")

-- Critical_points(f, a, b, opts) -> array of {x, f_x, kind}, status,
-- diagnostics
--
-- Every point in the CLOSED interval [a,b] where f'(x) = 0, each
-- classified as "MINIMUM", "MAXIMUM", or "NEITHER" (an inflection point
-- like x^3 at 0, where the derivative touches zero but does not change
-- sign) -- reusing critical_points.lua's own scan-then-classify pipeline
-- directly, the exact same one Global_minimum/Global_maximum in
-- optimization_numeric.lua already use to build their own candidate
-- list.
--
-- Same honest limitation as that scan: a coarse enough sample grid
-- (opts.samples, default 50) can still miss a very narrow feature
-- entirely between two samples -- a sampling method, not a certified
-- one.
--
-- Does NOT detect a critical point where f itself is not differentiable
-- (e.g. |x| at 0) -- this only ever looks at where f' EXISTS and equals
-- zero; a kink where f' doesn't exist at all is a different, harder
-- detection problem this function does not attempt.
local function Critical_points(f, a, b, opts)
    opts = opts or {}
    assert(type(opts) == "table", "Critical_points: opts must be a table.")
    f = CompileExpression.Resolve_function(f, "Critical_points: f")
    Core.assert_finite_number(a, "Critical_points: a")
    Core.assert_finite_number(b, "Critical_points: b")
    assert(a < b, "Critical_points: a must be less than b.")
    local samples = opts.samples or 50
    Core.assert_positive_integer(samples, "Critical_points: opts.samples")

    local g = CriticalPoints.derivative_fn(f)
    local xs = CriticalPoints.scan_critical_points(f, a, b, samples)
    table.sort(xs)

    local points = {}
    for _, x in ipairs(xs) do
        points[#points + 1] = {x = x, f_x = f(x), kind = CriticalPoints.classify_at(g, x)}
    end
    return points, "SUCCESS", {samples = samples}
end

-- Claude: Local_maxima/Local_minima are thin filters over Critical_points
-- -- ONE scan, not two independently-maintained searches that could
-- silently drift apart.
local function filter_by_kind(f, a, b, opts, wanted_kind)
    local points, status, diag = Critical_points(f, a, b, opts)
    local filtered = {}
    for _, p in ipairs(points) do
        if p.kind == wanted_kind then
            filtered[#filtered + 1] = p
        end
    end
    return filtered, status, diag
end

-- Local_maxima(f, a, b, opts) -> array of {x, f_x, kind="MAXIMUM"},
-- status, diagnostics
local function Local_maxima(f, a, b, opts)
    return filter_by_kind(f, a, b, opts, "MAXIMUM")
end

-- Local_minima(f, a, b, opts) -> array of {x, f_x, kind="MINIMUM"},
-- status, diagnostics
local function Local_minima(f, a, b, opts)
    return filter_by_kind(f, a, b, opts, "MINIMUM")
end

-- Claude: partitions [a,b] at every critical point found (the only
-- places f' can change sign for a well-behaved f), then samples f' at
-- each sub-interval's OWN midpoint to decide that whole sub-interval's
-- sign -- valid precisely because critical points are where sign changes
-- can happen, so between two consecutive ones the sign cannot flip.
local function monotonicity_intervals(f, a, b, opts, want_increasing)
    local points = Critical_points(f, a, b, opts)
    local boundaries = {a}
    for _, p in ipairs(points) do
        boundaries[#boundaries + 1] = p.x
    end
    boundaries[#boundaries + 1] = b

    local g = CriticalPoints.derivative_fn(f)
    local intervals = {}
    for i = 1, #boundaries - 1 do
        local left, right = boundaries[i], boundaries[i + 1]
        if right > left then
            local mid_slope = g((left + right) / 2)
            if (want_increasing and mid_slope > 0) or (not want_increasing and mid_slope < 0) then
                intervals[#intervals + 1] = {left, right}
            end
        end
    end
    return intervals, "SUCCESS", {samples = (opts or {}).samples or 50}
end

-- Increasing_intervals(f, a, b, opts) -> array of {left, right} pairs,
-- status, diagnostics
--
-- Every maximal sub-interval of [a,b] on which f' > 0 (strictly
-- increasing), found by partitioning at the critical points from
-- Critical_points and sampling f' at each piece's midpoint.
local function Increasing_intervals(f, a, b, opts)
    opts = opts or {}
    assert(type(opts) == "table", "Increasing_intervals: opts must be a table.")
    f = CompileExpression.Resolve_function(f, "Increasing_intervals: f")
    Core.assert_finite_number(a, "Increasing_intervals: a")
    Core.assert_finite_number(b, "Increasing_intervals: b")
    assert(a < b, "Increasing_intervals: a must be less than b.")
    return monotonicity_intervals(f, a, b, opts, true)
end

-- Decreasing_intervals(f, a, b, opts) -> array of {left, right} pairs,
-- status, diagnostics
--
-- The mirror image of Increasing_intervals -- every maximal sub-interval
-- of [a,b] on which f' < 0.
local function Decreasing_intervals(f, a, b, opts)
    opts = opts or {}
    assert(type(opts) == "table", "Decreasing_intervals: opts must be a table.")
    f = CompileExpression.Resolve_function(f, "Decreasing_intervals: f")
    Core.assert_finite_number(a, "Decreasing_intervals: a")
    Core.assert_finite_number(b, "Decreasing_intervals: b")
    assert(a < b, "Decreasing_intervals: a must be less than b.")
    return monotonicity_intervals(f, a, b, opts, false)
end

-- Tangent_line(f, x0) -> slope, intercept, L
--
-- The line tangent to f at x0: y = slope*x + intercept, where
-- slope = f'(x0) (exact, via Derivative_at_point) and
-- intercept = f(x0) - slope*x0. L is a plain callable Lua function,
-- L(x) = slope*x + intercept -- the same thing as the LINEARIZATION of f
-- at x0 (they are mathematically identical; this is deliberately not a
-- separate "Linearization" function that would just duplicate this one).
local function Tangent_line(f, x0)
    f = CompileExpression.Resolve_function(f, "Tangent_line: f")
    Core.assert_finite_number(x0, "Tangent_line: x0")
    local slope, value = DifNumeric.Derivative_at_point(f, x0)
    local intercept = value - slope * x0
    local L = function(x) return slope * x + intercept end
    return slope, intercept, L
end

-- Differential_approximation(f, x0, delta_x) -> approx, exact, error
--
-- The first-order (linear/differential) approximation of
-- f(x0 + delta_x), using only f(x0) and f'(x0): approx = f(x0) +
-- f'(x0)*delta_x -- the classical "df ~ f'(x0)*dx" approximation.
-- Returns the EXACT f(x0+delta_x) and the signed error too (approx -
-- exact), so a caller can see directly how the approximation degrades as
-- delta_x grows, rather than trusting the approximation blindly.
local function Differential_approximation(f, x0, delta_x)
    f = CompileExpression.Resolve_function(f, "Differential_approximation: f")
    Core.assert_finite_number(x0, "Differential_approximation: x0")
    Core.assert_finite_number(delta_x, "Differential_approximation: delta_x")
    local slope, value = DifNumeric.Derivative_at_point(f, x0)
    local approx = value + slope * delta_x
    local exact = f(x0 + delta_x)
    return approx, exact, approx - exact
end

return {
    Critical_points = Critical_points,
    Local_maxima = Local_maxima,
    Local_minima = Local_minima,
    Increasing_intervals = Increasing_intervals,
    Decreasing_intervals = Decreasing_intervals,
    Tangent_line = Tangent_line,
    Differential_approximation = Differential_approximation,
}
