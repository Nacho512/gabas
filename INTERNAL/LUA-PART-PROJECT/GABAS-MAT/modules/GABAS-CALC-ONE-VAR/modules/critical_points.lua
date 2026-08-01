-- [Critical Points]
-- Internal plumbing shared by optimization_numeric.lua and
-- derivative_analysis.lua: finding where f'(x) = 0 and classifying what
-- kind of critical point each one is. Not part of the public API (same
-- status as dual.lua, quadrature.lua, tanh_sinh.lua -- excluded from
-- GABAS_CALC_ONE_VAR.lua's aggregator list).
--
-- Claude: split out of optimization_numeric.lua (where this used to live
-- as private local functions) once derivative_analysis.lua needed the
-- SAME capability -- finding critical points is more primitive than
-- "which one optimizes f" (optimization compares the critical points
-- this module finds; it doesn't need to know how they were found), so
-- the dependency direction is optimization_numeric.lua depending on
-- THIS file, not the reverse, and not two independently-maintained
-- copies of the same scanning logic in two files.
local Core = require("GABAS-CALC-ONE-VAR.modules.core")
local DifNumeric = require("GABAS-CALC-ONE-VAR.modules.dif_numeric")
local RootsNumeric = require("GABAS-CALC-ONE-VAR.modules.roots_numeric")

-- Claude: g(x) = f'(x), computed exactly via Derivative_at_point (Dual
-- numbers) rather than a finite difference. f must therefore be
-- Dual-compatible (a symbolic-expression string, auto-compiled, or an
-- already-compiled function built only from Compile_expression's own
-- Dual-aware primitives) -- same requirement Newton in roots_numeric.lua
-- already has, for the same reason.
local function derivative_fn(f)
    return function(x) return (DifNumeric.Derivative_at_point(f, x)) end
end

-- Claude: classifies a KNOWN critical point x (already g(x) approx 0) by
-- checking the sign of g = f' immediately to either side -- MINIMUM
-- (- to +), MAXIMUM (+ to -), or NEITHER (an inflection point like x^3
-- at 0, where f'=0 but does not change sign). The probe distance uses
-- the same cube-root-of-machine-epsilon scaling as
-- Verify_derivative_at_point in dif_numeric.lua -- small enough to probe
-- genuinely local behavior, large enough to stay clear of floating-point
-- noise. Split out as its own function (rather than folded into
-- find_critical_point_and_classify below) so a caller that already has a
-- precise x (e.g. derivative_analysis.lua's Critical_points, working
-- from scan_critical_points' own already-refined results) can classify
-- it directly, without re-running a root search that would just
-- converge back to the same point.
local function classify_at(g, x)
    local eps = (Core.MACHINE_EPSILON ^ (1 / 3)) * math.max(1, math.abs(x))
    local slope_left = g(x - eps)
    local slope_right = g(x + eps)
    if slope_left < 0 and slope_right > 0 then
        return "MINIMUM"
    elseif slope_left > 0 and slope_right < 0 then
        return "MAXIMUM"
    end
    return "NEITHER"
end

-- Claude: finds ONE critical point near x0 (root-finding on f' via
-- Secant, seeded with x0 and a nearby second point) and classifies it
-- via classify_at above.
local function find_critical_point_and_classify(f, x0, opts)
    local g = derivative_fn(f)
    local step = opts.initial_step or 0.01 * math.max(1, math.abs(x0))
    local x_star, _, status, diag = RootsNumeric.Secant(g, x0, x0 + step, opts)
    if status ~= "SUCCESS" then
        return x_star, nil, nil, status, diag
    end
    return x_star, f(x_star), classify_at(g, x_star), "SUCCESS", diag
end

-- Claude: samples f' at `samples` evenly spaced points across [a,b] and
-- refines every detected sign change (via Secant on f') into a precise
-- critical x. Deliberately does NOT classify each one here -- callers
-- that need the kind (derivative_analysis.lua's Critical_points) run
-- find_critical_point_and_classify on each result; callers that only
-- need the raw x values to compare f-values (optimization_numeric.lua's
-- Global_minimum/maximum) skip that extra work. A sample that lands (to
-- floating-point precision) exactly on a critical point is also caught
-- directly, rather than relying on a sign change either side of it.
--
-- Honest limitation, stated plainly rather than hidden: a coarse enough
-- sample grid can still miss a very narrow feature entirely between two
-- samples -- raising `samples` improves the resolution at the cost of
-- more evaluations, but this is a sampling method, not a certified one
-- (that would need interval arithmetic, a different, future technique).
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

return {
    derivative_fn = derivative_fn,
    classify_at = classify_at,
    find_critical_point_and_classify = find_critical_point_and_classify,
    scan_critical_points = scan_critical_points,
}
