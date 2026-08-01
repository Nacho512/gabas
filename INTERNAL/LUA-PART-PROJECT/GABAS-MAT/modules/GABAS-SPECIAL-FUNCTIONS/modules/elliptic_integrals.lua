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

-- Claude: K(m), E(m), F(phi,m), and E(phi,m) below are all standard
-- Carlson-form reductions (DLMF 19.25.1, 19.25.5, 19.25.9) -- copied
-- from that reference, not derived from scratch, same discipline as
-- the RC/RF/RJ pseudocode above. They fall directly out of the third-
-- kind reduction the reference document derives (Appendix A): F is
-- exactly Pi's own RF term alone (Pi(0;phi|m) = F(phi|m), the n=0
-- special case the document itself dispatches to), and E uses the
-- identity RD(x,y,z) = RJ(x,y,z,z) the document's own test plan lists
-- (M8) -- no separate RD kernel is implemented; RJ(x,y,1,1) already IS
-- RD(x,y,1) by that identity.
--
--   F(phi|m)  = s * RF(c^2, y, 1)
--   E(phi|m)  = s * RF(c^2, y, 1) - (m*s^3/3) * RJ(c^2, y, 1, 1)
--   K(m)      = F(pi/2|m) = RF(0, 1-m, 1)
--   E(m)      = E(pi/2|m) = RF(0, 1-m, 1) - (m/3) * RJ(0, 1-m, 1, 1)
--
-- where s = sin(phi), c = cos(phi), y = 1 - m*s^2.
--
-- Domain: m < 1, |phi| <= pi/2 for the incomplete forms (m < 1 alone
-- for the complete forms, phi = pi/2 fixed). Deliberately conservative,
-- same "verified, not theoretical" discipline as every other domain
-- boundary in this project -- and not an arbitrary restriction: for
-- m < 1, y = 1 - m*sin^2(phi) >= 1 - m > 0 STRICTLY for every phi in
-- this range, so c^2 is the only argument that can ever be zero (only
-- at the endpoints +-pi/2), never y too -- Carlson_rf's "at most one
-- zero" precondition can never be violated on this domain. m >= 1
-- introduces real additional cases (K diverges as m approaches 1;
-- m > 1 needs its own radicand-root bookkeeping, exactly the kind of
-- fringe-case machinery the reference document spends a full chapter
-- on for the harder third-kind integral) -- deferred, real future
-- work, not attempted here.
--
-- Verified directly against mpmath.ellipf/ellipe/ellipk (independent,
-- arbitrary-precision references): relative error stays within ~1e-15
-- (essentially machine epsilon, inherited directly from the RF/RJ
-- kernels' own verified accuracy) across a wide sweep of phi and m,
-- including m very negative.
local function Elliptic_f(phi, m)
    Core.assert_finite_number(phi, "Elliptic_f: phi")
    Core.assert_finite_number(m, "Elliptic_f: m")
    assert(phi >= -math.pi / 2 and phi <= math.pi / 2, "Elliptic_f: phi must be in [-pi/2, pi/2].")
    assert(m < 1, "Elliptic_f: m must be < 1.")
    local s, c = math.sin(phi), math.cos(phi)
    local y = 1 - m * s * s
    return s * Carlson_rf(c * c, y, 1)
end

local function Elliptic_e_incomplete(phi, m)
    Core.assert_finite_number(phi, "Elliptic_e_incomplete: phi")
    Core.assert_finite_number(m, "Elliptic_e_incomplete: m")
    assert(phi >= -math.pi / 2 and phi <= math.pi / 2, "Elliptic_e_incomplete: phi must be in [-pi/2, pi/2].")
    assert(m < 1, "Elliptic_e_incomplete: m must be < 1.")
    local s, c = math.sin(phi), math.cos(phi)
    local y = 1 - m * s * s
    return s * Carlson_rf(c * c, y, 1) - (m * s ^ 3 / 3) * Carlson_rj(c * c, y, 1, 1)
end

local function Elliptic_k(m)
    Core.assert_finite_number(m, "Elliptic_k: m")
    assert(m < 1, "Elliptic_k: m must be < 1.")
    return Carlson_rf(0, 1 - m, 1)
end

local function Elliptic_e(m)
    Core.assert_finite_number(m, "Elliptic_e: m")
    assert(m < 1, "Elliptic_e: m must be < 1.")
    return Carlson_rf(0, 1 - m, 1) - (m / 3) * Carlson_rj(0, 1 - m, 1, 1)
end

-- Elliptic_pi(n, phi, m) -> Pi(n;phi|m), the incomplete elliptic
-- integral of the third kind:
--
--   Pi(n;phi|m) = integral from 0 to phi of
--                 dtheta / [(1 - n*sin^2(theta)) * sqrt(1 - m*sin^2(theta))]
--
-- Claude: the function the reference document Nacho supplied is
-- entirely about -- read in full before any of this file was written.
-- Its own recommended production design (Section "Executive
-- recommendation", Chapter 3) is exactly this Carlson reduction:
--
--   Pi(n;phi|m) = s*RF(c^2, y, 1) + (n*s^3/3)*RJ(c^2, y, 1, p)
--
-- where s = sin(phi), c = cos(phi), y = 1-m*s^2, p = 1-n*s^2 -- the
-- SAME x,y,z Elliptic_f/Elliptic_e_incomplete already use, plus one
-- more RJ call with the third-kind denominator p as its fourth
-- argument. Derived from the substitution t=sin(theta) in Appendix A of
-- the reference document.
--
-- Domain: m < 1, n < 1, |phi| <= pi/2 -- deliberately restricted to
-- what the document itself calls the "globally nonsingular periodic
-- regime" (Section 1.4: "If m<1 and n<1, there are no real
-- singularities in a period and the complete integral is finite").
-- Not an arbitrary restriction: n < 1 guarantees
-- p = 1 - n*sin^2(theta) > 0 STRICTLY for every theta in this range
-- (mirroring exactly how m < 1 already keeps y > 0 strictly for
-- Elliptic_f/Elliptic_e_incomplete above), so Carlson_rj's p > 0
-- precondition can never be violated here, and there is no pole to
-- cross, classify, or take a principal value around. The document
-- spends most of its length (domain classification, POLE_CROSSED /
-- NONREAL_PATH statuses, Cauchy principal values, periodic amplitude
-- reduction past +-pi/2, Bulirsch/Fukushima alternatives) on exactly
-- the n >= 1 / m >= 1 / huge-|phi| cases this domain excludes --
-- deferred, real future work, not attempted here, same discipline as
-- every other domain boundary in this project.
--
-- No special-case dispatch (n=0 -> F, m=0 -> elementary, n=m -> E) is
-- needed for correctness on this domain, unlike what the document's
-- general-domain dispatcher requires -- checked directly, not assumed:
-- at n=0 the RJ term is multiplied by n=0 and vanishes on its own,
-- exactly reproducing Elliptic_f(phi,m); m=0 and n=m both evaluate
-- correctly straight through the general Carlson formula, with no
-- cancellation or singularity in this restricted domain (the document
-- itself notes the n=m closed form is only ever a cross-check
-- alternative, not a numerical necessity, away from m=1).
--
-- Verified directly against mpmath.ellippi (an independent,
-- arbitrary-precision reference): relative error at the few-1e-16
-- (machine epsilon) level across 300 random (n,phi,m) points in this
-- domain, and against 16 of the reference document's own 20
-- high-precision test vectors (Table 14.1) that fall inside it --
-- including T16/T17/T18, which stress m and n individually approaching
-- 1 to within 1e-9/1e-12. (T06, T07, T10, and T14 fall outside this
-- domain -- n or m >= 1, or |phi| needing periodic reduction -- and are
-- correctly rejected rather than silently mishandled.)
local function Elliptic_pi(n, phi, m)
    Core.assert_finite_number(n, "Elliptic_pi: n")
    Core.assert_finite_number(phi, "Elliptic_pi: phi")
    Core.assert_finite_number(m, "Elliptic_pi: m")
    assert(phi >= -math.pi / 2 and phi <= math.pi / 2, "Elliptic_pi: phi must be in [-pi/2, pi/2].")
    assert(m < 1, "Elliptic_pi: m must be < 1.")
    assert(n < 1, "Elliptic_pi: n must be < 1 (n >= 1 introduces third-kind pole-crossing cases -- a deferred future increment).")
    local s, c = math.sin(phi), math.cos(phi)
    local y = 1 - m * s * s
    local p = 1 - n * s * s
    return s * Carlson_rf(c * c, y, 1) + (n * s ^ 3 / 3) * Carlson_rj(c * c, y, 1, p)
end

return {
    Carlson_rc = Carlson_rc,
    Carlson_rf = Carlson_rf,
    Carlson_rj = Carlson_rj,
    Elliptic_f = Elliptic_f,
    Elliptic_e_incomplete = Elliptic_e_incomplete,
    Elliptic_k = Elliptic_k,
    Elliptic_e = Elliptic_e,
    Elliptic_pi = Elliptic_pi,
}
