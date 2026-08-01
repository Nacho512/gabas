-- [Tanh Sinh]
-- The double-exponential (DE) change of variable used by tanh-sinh
-- quadrature: u(t) = tanh(pi/2 * sinh(t)) maps the whole real line in t
-- onto the open interval (-1,1) in u, with a Jacobian weight that decays
-- DOUBLY exponentially as |t| grows -- fast enough to overwhelm an
-- algebraic (integrable) singularity of the integrand near u = +-1,
-- which is why this method tolerates endpoint singularities that an
-- ordinary rule (Gauss-Kronrod included) cannot.
--
-- Claude: this file holds only the raw node/weight mechanics for ONE
-- evaluation point -- De_map below, and the private real_sinh/real_cosh
-- it's built from. The orchestrating piece (level-doubling trapezoidal
-- sums, convergence checking) lives in integration_numeric.lua's
-- Integral_tanh_sinh, mirroring how quadrature.lua (one Gauss-Kronrod
-- panel) is kept separate from Integral_definite (the adaptive loop
-- built on top of it).
--
-- Real_sinh/real_cosh are defined locally rather than relying on Lua's
-- own math.sinh/cosh -- those were made optional (LUA_COMPAT_5_2) as of
-- Lua 5.3 and aren't guaranteed to exist in every build, including
-- whatever Lua LuaTeX embeds, the actual eventual target for this
-- module (same reasoning already applied in compile_expression.lua and
-- GABAS-LINAL's complex.lua).
local function real_sinh(t) return (math.exp(t) - math.exp(-t)) / 2 end
local function real_cosh(t) return (math.exp(t) + math.exp(-t)) / 2 end

local HALF_PI = math.pi / 2

-- De_map(t) -> delta, w
--
-- delta = 1 - tanh(|pi/2 * sinh(t)|), always in (0, 2] -- the distance
-- from the endpoint being approached (small, for large |t|).
-- w = du/dt = (pi/2 * cosh(t)) / cosh(pi/2*sinh(t))^2, the Jacobian
-- weight for the trapezoidal sum in t.
--
-- Claude: delta is computed DIRECTLY as 2/(e^(2*arg)+1), NEVER by first
-- forming u = tanh(arg) and then subtracting 1-u -- an earlier version
-- of this file did exactly that, and it was a real, confirmed bug.
-- (e2x-1)/(e2x+1) computes the DIFFERENCE of two huge, nearly-equal
-- numbers, which loses the tiny distinguishing information to
-- floating-point rounding once e2x exceeds about 1e15 -- saturating u to
-- EXACTLY +-1 around arg ~ 17, even though the true value there is still
-- meaningfully less than 1 (1 - ~1.5e-15) and the true x it corresponds
-- to is still meaningfully inside (a,b), not AT the endpoint. Since a
-- node this close to a genuine endpoint singularity can still carry a
-- non-negligible contribution (traced directly: for integral of x^-0.5
-- on [0,1], the weight*integrand product at the last node before this
-- bug's premature cutoff was still ~9e-6, not remotely negligible), that
-- premature saturation produced a hard, refinement-proof ~1e-8 error
-- that no amount of step-size halving could fix (the SAME tail was
-- always missing, at every level). Computing delta = 2/(e2arg+1)
-- directly is an ordinary division, not a cancellation -- it stays
-- accurate all the way out to where e2arg itself overflows (arg ~ 355,
-- roughly 20x further than the old u-based cutoff), which is also
-- comfortably before cosh(arg) (needed for w) would overflow.
--
-- The caller combines delta with the known endpoints: for t>0
-- (approaching b), x = b - half_range*delta; for t<0 (approaching a),
-- x = a + half_range*delta. delta depends only on |t|, so ONE call to
-- De_map(|t|) serves both the +t and -t nodes of a level.
--
-- Returns nil, nil -- meaning "this node contributes nothing, do not
-- evaluate the integrand here" -- once delta has genuinely underflowed
-- to 0 (no floating-point value can represent the node any closer to the
-- endpoint -- technical document dangers #2/#3) or the weight has
-- (danger #4).
local function De_map(t)
    local abs_t = math.abs(t)
    local sinh_t = real_sinh(abs_t)
    local cosh_t = real_cosh(abs_t)
    local arg = HALF_PI * sinh_t
    local e2arg = math.exp(2 * arg)
    local delta
    if e2arg == math.huge then
        delta = 0
    else
        delta = 2 / (e2arg + 1)
    end
    if delta <= 0 then
        return nil, nil
    end
    local cosh_arg = real_cosh(arg)
    local w = (HALF_PI * cosh_t) / (cosh_arg * cosh_arg)
    if w ~= w or w == 0 then
        return nil, nil
    end
    return delta, w
end

return {
    De_map = De_map,
}
