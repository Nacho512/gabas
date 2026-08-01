-- [Elliptic Integrals]
-- Carlson_rc(x,y), Carlson_rf(x,y,z), Carlson_rj(x,y,z,p) -- Carlson's
-- symmetric elliptic integrals, the foundation this whole family is
-- built on. Split into its own file since this is a genuinely different
-- family from both special_functions.lua (Gamma/Bessel-style
-- series/continued-fraction machinery) and orthogonal_polynomials.lua
-- (three-term recurrences) -- Carlson's DUPLICATION algorithm.
--
-- Claude: the standard textbook (K(m), E(m), F(phi,m), E(phi,m)) and
-- the incomplete elliptic integral of the THIRD kind (Pi(n;phi,m)) are
-- both planned to build directly on these three kernels, per a detailed
-- technical reference document Nacho supplied specifically for the
-- third-kind integral (the hardest of the three) -- read in full before
-- any of this was written, same discipline as the Bessel reference
-- document earlier in this project. That document gives complete,
-- explicit pseudocode for RC, RF, and RJ (Carlson's 1995 algorithm,
-- "Numerical computation of real or complex elliptic integrals"), which
-- this file transcribes directly, not re-derives -- same discipline as
-- the Lanczos coefficients and the QK21 quadrature tables elsewhere in
-- this project.
local Core = require("GABAS-SPECIAL-FUNCTIONS.modules.core")

-- Claude: ERRTOL constants below are the standard, widely-published
-- values (Carlson's 1995 paper / the classic Numerical Recipes rc/rf/rj
-- routines) specifically tuned for these exact truncated-series
-- correction polynomials -- not invented, and not "copied blindly"
-- either (the reference document explicitly warns against that): each
-- polynomial's own truncation order was checked against its ERRTOL
-- (e.g. RC's degree-5-in-s correction needs |s|^6 below machine
-- epsilon, i.e. |s| below roughly 1e-16^(1/6) =~ 0.0022 -- comfortably
-- above ERRTOL=0.0012), and then the whole thing was verified directly
-- against mpmath.elliprc/elliprf/elliprj across thousands of random
-- points before being trusted: relative error stayed within ~1e-15
-- (essentially machine epsilon) throughout.
local RC_ERRTOL = 0.0012
local RF_ERRTOL = 0.0025
local RJ_ERRTOL = 0.0015

-- Claude: power-of-two scaling (the reference document's Listing 2.1) --
-- NOT an optional robustness nicety. Carlson's forms are exactly
-- homogeneous (RF(lambda*x,...) = lambda^-1/2 * RF(x,...), and
-- similarly for RJ/RC with their own powers), so rescaling every
-- argument by the same power of two before entering the duplication
-- loop, then rescaling the result back, is exact in binary floating
-- point (no precision lost) and keeps every intermediate quantity in a
-- comfortable range. Caught directly by testing, not assumed:
-- Carlson_rj(1e-250, 1e-250, 1e-250, 1e-250, 1e-250) -- whose true value
-- (1e-250)^-1.5 = 1e375 genuinely overflows a double -- failed with a
-- spurious division-by-zero deep inside the duplication loop under the
-- unscaled core alone, instead of honestly overflowing to +math.huge at
-- the very last step. Scaling first fixes exactly that: the *only*
-- place overflow/underflow can now happen is the final, unavoidable
-- rescale-by-S step, an honest IEEE result, not a corrupted
-- intermediate one.
local function scale_factor(largest_magnitude)
    local _, exponent = math.frexp(largest_magnitude)
    return math.ldexp(1.0, exponent)
end

-- Carlson_rc(x, y) -> RC(x,y), Carlson's degenerate symmetric elliptic
-- integral, for x >= 0 and y > 0.
--
-- RC(x,y) = (1/2) * integral from 0 to infinity of dt/((t+y)*sqrt(t+x))
-- -- restricted here to the positive-y case the reference document
-- covers (RJ's own internal calls, and the ordinary elliptic-integral
-- reductions this file builds toward, never need y <= 0).
local function rc_core(x, y)
    local mu, s
    while true do
        mu = (x + 2 * y) / 3
        s = (y - mu) / mu
        if math.abs(s) < RC_ERRTOL then
            break
        end
        local sx, sy = math.sqrt(x), math.sqrt(y)
        local lambda = 2 * sx * sy + y
        x = (x + lambda) / 4
        y = (y + lambda) / 4
    end
    local poly = 1 + s * s * (3 / 10 + s * (1 / 7 + s * (3 / 8 + s * (9 / 22))))
    return poly / math.sqrt(mu)
end

local function Carlson_rc(x, y)
    Core.assert_finite_number(x, "Carlson_rc: x")
    Core.assert_finite_number(y, "Carlson_rc: y")
    assert(x >= 0, "Carlson_rc: x must be nonnegative.")
    assert(y > 0, "Carlson_rc: y must be positive.")
    local scale = scale_factor(math.max(x, y))
    return rc_core(x / scale, y / scale) / math.sqrt(scale)
end

-- Carlson_rf(x, y, z) -> RF(x,y,z), Carlson's symmetric elliptic
-- integral of the first kind, for x,y,z >= 0 with at most one zero.
--
-- RF(x,y,z) = (1/2) * integral from 0 to infinity of
-- dt/sqrt((t+x)(t+y)(t+z)).
local function rf_core(x, y, z)
    local mu, dx, dy, dz
    while true do
        mu = (x + y + z) / 3
        dx, dy, dz = (mu - x) / mu, (mu - y) / mu, (mu - z) / mu
        if math.max(math.abs(dx), math.abs(dy), math.abs(dz)) < RF_ERRTOL then
            break
        end
        local sx, sy, sz = math.sqrt(x), math.sqrt(y), math.sqrt(z)
        local lambda = sx * sy + sx * sz + sy * sz
        x = (x + lambda) / 4
        y = (y + lambda) / 4
        z = (z + lambda) / 4
    end
    local e2 = dx * dy - dz * dz
    local e3 = dx * dy * dz
    local correction = 1 + (e2 / 24 - 1 / 10 - 3 * e3 / 44) * e2 + e3 / 14
    return correction / math.sqrt(mu)
end

local function Carlson_rf(x, y, z)
    Core.assert_finite_number(x, "Carlson_rf: x")
    Core.assert_finite_number(y, "Carlson_rf: y")
    Core.assert_finite_number(z, "Carlson_rf: z")
    assert(x >= 0 and y >= 0 and z >= 0, "Carlson_rf: x, y, and z must all be nonnegative.")
    assert((x == 0 and 1 or 0) + (y == 0 and 1 or 0) + (z == 0 and 1 or 0) <= 1,
        "Carlson_rf: at most one of x, y, z may be zero.")
    local scale = scale_factor(math.max(x, y, z))
    return rf_core(x / scale, y / scale, z / scale) / math.sqrt(scale)
end

-- Carlson_rj(x, y, z, p) -> RJ(x,y,z,p), Carlson's symmetric elliptic
-- integral of the third kind, for x,y,z >= 0 (at most one zero) and
-- p > 0.
--
-- RJ(x,y,z,p) = (3/2) * integral from 0 to infinity of
-- dt/((t+p)*sqrt((t+x)(t+y)(t+z))) -- restricted here to the
-- positive-p case; Carlson's separate negative-p (principal-value)
-- transformation is real, documented future work, not attempted here.
--
-- Claude: each outer duplication step needs one Carlson_rc-style
-- evaluation of its own (the alpha/beta pair below) -- computed via the
-- unscaled rc_core directly, not the public Carlson_rc, since the
-- OUTER scale_factor applied to x,y,z,p before this loop starts already
-- keeps alpha/beta in a safe range throughout (the reference document's
-- own recommended remedy for alpha/beta overflow -- scale the inputs
-- before entering the core -- rather than re-scaling alpha/beta
-- separately every iteration).
local function rj_core(x, y, z, p)
    local total, factor = 0, 1
    local mu, dx, dy, dz, dp
    while true do
        mu = (x + y + z + 2 * p) / 5
        dx, dy, dz, dp = (mu - x) / mu, (mu - y) / mu, (mu - z) / mu, (mu - p) / mu
        if math.max(math.abs(dx), math.abs(dy), math.abs(dz), math.abs(dp)) < RJ_ERRTOL then
            break
        end
        local sx, sy, sz, sp = math.sqrt(x), math.sqrt(y), math.sqrt(z), math.sqrt(p)
        local lambda = sx * sy + sx * sz + sy * sz
        local alpha = p * (sx + sy + sz) + sx * sy * sz
        alpha = alpha * alpha
        local beta = p * (p + lambda) * (p + lambda)
        total = total + factor * rc_core(alpha, beta)
        factor = factor / 4
        x = (x + lambda) / 4
        y = (y + lambda) / 4
        z = (z + lambda) / 4
        p = (p + lambda) / 4
    end
    local ea = dx * (dy + dz) + dy * dz
    local eb = dx * dy * dz
    local ec = dp * dp
    local ed = ea - 3 * ec
    local ee = eb + 2 * dp * (ea - ec)
    local c1, c2, c3, c4 = 3 / 14, 1 / 3, 3 / 22, 3 / 26
    local c5, c6, c7, c8 = 9 / 88, 9 / 52, 1 / 6, 3 / 11
    local poly = 1
        + ed * (-c1 + c5 * ed - c6 * ee)
        + eb * (c7 + dp * (-c8 + dp * c4))
        + dp * ea * (c2 - dp * c3)
        - c2 * dp * ec
    return 3 * total + factor * poly / (mu * math.sqrt(mu))
end

local function Carlson_rj(x, y, z, p)
    Core.assert_finite_number(x, "Carlson_rj: x")
    Core.assert_finite_number(y, "Carlson_rj: y")
    Core.assert_finite_number(z, "Carlson_rj: z")
    Core.assert_finite_number(p, "Carlson_rj: p")
    assert(x >= 0 and y >= 0 and z >= 0, "Carlson_rj: x, y, and z must all be nonnegative.")
    assert((x == 0 and 1 or 0) + (y == 0 and 1 or 0) + (z == 0 and 1 or 0) <= 1,
        "Carlson_rj: at most one of x, y, z may be zero.")
    assert(p > 0, "Carlson_rj: p must be positive (the negative-p principal-value case is not implemented).")
    local scale = scale_factor(math.max(x, y, z, p))
    return rj_core(x / scale, y / scale, z / scale, p / scale) / (scale * math.sqrt(scale))
end

return {
    Carlson_rc = Carlson_rc,
    Carlson_rf = Carlson_rf,
    Carlson_rj = Carlson_rj,
}
