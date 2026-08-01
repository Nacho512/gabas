-- [Roots Numeric]
-- Numerical root-finding for a REAL one-variable function: Bisection,
-- Newton-Raphson, and Secant. Named "roots_numeric" (not just "roots") to
-- leave the name "roots_symbolic" free for a future closed-form/symbolic
-- root solver (quadratic formula, polynomial factoring, ...) -- that is
-- CAS-territory, a different, later phase of this project, exactly like
-- the differentiation_symbolic/integration_symbolic split already visible
-- in the reference architecture document this project consulted (and
-- deliberately did not adopt wholesale -- see the project's own
-- discussion of that document).
--
-- Claude: every function here requires f to be REAL-valued at every point
-- it is evaluated -- f(x) must be a plain Lua number, not a Complex
-- value. This is not an incidental restriction: bisection's sign test and
-- secant's ordering both presuppose a real-valued, real-ordered f, and
-- Core.assert_finite_number below rejects a Complex result outright
-- (type(x) == "number" fails for a Complex table) rather than silently
-- misbehaving on an ordering comparison that isn't mathematically
-- meaningful for complex numbers. Consistent with this whole package's
-- scope: REAL one-variable calculus only.
local Core = require("GABAS-CALC-ONE-VAR.modules.core")
local CompileExpression = require("GABAS-CALC-ONE-VAR.modules.compile_expression")
local DifNumeric = require("GABAS-CALC-ONE-VAR.modules.dif_numeric")

-- Bisection(f, a, b, opts) -> root, estimated_error, status, diagnostics
--
-- The classical bracketing method: f(a) and f(b) must have opposite
-- signs (a genuine sign change on [a,b], guaranteed by the intermediate
-- value theorem to contain at least one root for a continuous f).
-- Halves the bracket every iteration, keeping whichever half still
-- contains the sign change. Always converges given a valid bracket --
-- the tradeoff for that guarantee is linear (not superlinear)
-- convergence.
--
-- opts (all optional): tol (default 1e-12, the bracket half-width at
-- which to stop), max_iterations (default 100).
--
-- estimated_error is the final bracket half-width -- a genuine, provable
-- bound on |root - true root| (not just a heuristic), unlike Newton's or
-- Secant's step-size-based estimate below.
local function Bisection(f, a, b, opts)
    opts = opts or {}
    assert(type(opts) == "table", "Bisection: opts must be a table.")
    f = CompileExpression.Resolve_function(f, "Bisection: f")
    Core.assert_finite_number(a, "Bisection: a")
    Core.assert_finite_number(b, "Bisection: b")
    assert(a ~= b, "Bisection: a and b must be distinct.")
    if b < a then a, b = b, a end

    local tol = opts.tol or 1e-12
    Core.assert_positive_number(tol, "Bisection: opts.tol")
    local max_iterations = opts.max_iterations or 100
    Core.assert_positive_integer(max_iterations, "Bisection: opts.max_iterations")

    local fa = f(a)
    Core.assert_finite_number(fa, "Bisection: f(a)")
    if fa == 0 then
        return a, 0, "SUCCESS", {iterations = 0, residual = 0}
    end
    local fb = f(b)
    Core.assert_finite_number(fb, "Bisection: f(b)")
    if fb == 0 then
        return b, 0, "SUCCESS", {iterations = 0, residual = 0}
    end
    assert((fa < 0) ~= (fb < 0),
        "Bisection: f(a) and f(b) must have opposite signs (no sign change detected on [a,b]).")

    local mid, fm
    for iteration = 1, max_iterations do
        mid = a + (b - a) / 2
        fm = f(mid)
        Core.assert_finite_number(fm, "Bisection: f(mid)")
        local half_width = (b - a) / 2
        if fm == 0 or half_width <= tol then
            return mid, half_width, "SUCCESS", {iterations = iteration, residual = math.abs(fm)}
        end
        if (fm < 0) == (fa < 0) then
            a, fa = mid, fm
        else
            b, fb = mid, fm
        end
    end

    return mid, (b - a) / 2, "MAX_ITERATIONS_REACHED", {iterations = max_iterations, residual = math.abs(fm)}
end

-- Newton(f, x0, opts) -> root, estimated_error, status, diagnostics
--
-- Newton-Raphson iteration: x_(n+1) = x_n - f(x_n)/f'(x_n). f'(x_n) is
-- NOT supplied by the caller and never approximated by finite
-- differences -- it comes from this package's own Derivative_at_point
-- (forward-mode automatic differentiation via Dual numbers), exact up to
-- floating-point rounding. This is a direct, deliberate reuse of already-
-- built, already-tested infrastructure: Newton's method and automatic
-- differentiation were built as two separate, independently verified
-- pieces, and combine here with no new derivative logic of their own.
--
-- opts (all optional): tol (default 1e-12, the step-size convergence
-- threshold), max_iterations (default 50), min_derivative (default
-- 1e-14, guards against dividing by a derivative that has collapsed to
-- (numerically) zero).
--
-- Converges quadratically near a SIMPLE root of a smooth f, but is not
-- guaranteed to converge at all from an arbitrary starting guess -- no
-- bracket is required, and none is checked. A residual already at (or
-- below) tol at the current iterate returns immediately, before ever
-- computing a step -- this is what makes x0 already being an exact root
-- (where f'(x0) might ALSO be 0, e.g. f(x)=x^2 at x0=0) succeed instead
-- of falsely tripping the min_derivative guard below.
local function Newton(f, x0, opts)
    opts = opts or {}
    assert(type(opts) == "table", "Newton: opts must be a table.")
    f = CompileExpression.Resolve_function(f, "Newton: f")
    Core.assert_finite_number(x0, "Newton: x0")

    local tol = opts.tol or 1e-12
    Core.assert_positive_number(tol, "Newton: opts.tol")
    local max_iterations = opts.max_iterations or 50
    Core.assert_positive_integer(max_iterations, "Newton: opts.max_iterations")
    local min_derivative = opts.min_derivative or 1e-14
    Core.assert_positive_number(min_derivative, "Newton: opts.min_derivative")

    local x = x0
    for iteration = 1, max_iterations do
        local derivative, value = DifNumeric.Derivative_at_point(f, x)
        Core.assert_finite_number(value, "Newton: f(x)")
        Core.assert_finite_number(derivative, "Newton: f'(x)")
        if math.abs(value) <= tol then
            return x, 0, "SUCCESS", {iterations = iteration - 1, residual = math.abs(value)}
        end
        if math.abs(derivative) < min_derivative then
            return x, nil, "DERIVATIVE_TOO_SMALL", {iterations = iteration - 1, residual = math.abs(value)}
        end
        local step = value / derivative
        local x_next = x - step
        if not Core.is_finite_number(x_next) then
            return x, nil, "DIVERGED", {iterations = iteration, residual = math.abs(value)}
        end
        if math.abs(step) <= tol then
            return x_next, math.abs(step), "SUCCESS", {iterations = iteration, residual = math.abs(value)}
        end
        x = x_next
    end

    return x, nil, "MAX_ITERATIONS_REACHED", {iterations = max_iterations}
end

-- Secant(f, x0, x1, opts) -> root, estimated_error, status, diagnostics
--
-- Approximates Newton's method's derivative with a finite difference
-- between the two most recent iterates, so no derivative (automatic or
-- otherwise) is needed -- useful when f is expensive to seed with a Dual
-- (or the caller would rather not depend on Derivative_at_point at all).
-- Superlinear convergence (order ~1.618, the golden ratio) near a simple
-- root -- slower than Newton's quadratic convergence, faster than
-- Bisection's linear convergence, and (like Newton) not guaranteed to
-- converge from an arbitrary pair of starting points.
--
-- opts (all optional): tol (default 1e-12), max_iterations (default 100).
local function Secant(f, x0, x1, opts)
    opts = opts or {}
    assert(type(opts) == "table", "Secant: opts must be a table.")
    f = CompileExpression.Resolve_function(f, "Secant: f")
    Core.assert_finite_number(x0, "Secant: x0")
    Core.assert_finite_number(x1, "Secant: x1")
    assert(x0 ~= x1, "Secant: x0 and x1 must be distinct.")

    local tol = opts.tol or 1e-12
    Core.assert_positive_number(tol, "Secant: opts.tol")
    local max_iterations = opts.max_iterations or 100
    Core.assert_positive_integer(max_iterations, "Secant: opts.max_iterations")

    local f0 = f(x0)
    Core.assert_finite_number(f0, "Secant: f(x0)")
    local f1 = f(x1)
    Core.assert_finite_number(f1, "Secant: f(x1)")

    for iteration = 1, max_iterations do
        if math.abs(f1) <= tol then
            return x1, 0, "SUCCESS", {iterations = iteration - 1, residual = math.abs(f1)}
        end
        if f1 == f0 then
            return x1, nil, "STALLED", {iterations = iteration - 1, residual = math.abs(f1)}
        end
        local x2 = x1 - f1 * (x1 - x0) / (f1 - f0)
        if not Core.is_finite_number(x2) then
            return x1, nil, "DIVERGED", {iterations = iteration, residual = math.abs(f1)}
        end
        local step = x2 - x1
        x0, f0 = x1, f1
        x1 = x2
        f1 = f(x1)
        Core.assert_finite_number(f1, "Secant: f(x)")
        if math.abs(step) <= tol then
            return x1, math.abs(step), "SUCCESS", {iterations = iteration, residual = math.abs(f1)}
        end
    end

    return x1, nil, "MAX_ITERATIONS_REACHED", {iterations = max_iterations, residual = math.abs(f1)}
end

return {
    Bisection = Bisection,
    Newton = Newton,
    Secant = Secant,
}
