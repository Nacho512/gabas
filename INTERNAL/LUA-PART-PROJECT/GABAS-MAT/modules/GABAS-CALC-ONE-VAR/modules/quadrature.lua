-- [Quadrature]
-- The 21-point Gauss-Kronrod rule (embedding a 10-point Gauss rule) on a
-- single subinterval -- the local building block Integral_definite (in
-- GABAS_CALC_ONE_VAR.lua) subdivides adaptively. This is the same rule
-- QUADPACK's QK21 uses, and the node/weight tables below are copied from
-- that reference (not re-derived) -- the technical specification this
-- module follows explicitly warns against improvising these constants
-- (section 6.4): "The constants should be taken from a validated
-- reference implementation rather than improvised."
--
-- Claude: sanity-checked before ANY use here -- sum(WGK) (each of the 10
-- paired entries counted twice, plus the center once) equals exactly 2,
-- and sum(WG) (each of the 5 paired entries counted twice) equals exactly
-- 2 as well, both to full double precision -- the two independent checks
-- a table transcribed correctly must satisfy (integrating the constant
-- function 1 over [-1,1] must give the length of the interval, 2, under
-- either rule). Also verified against known exact integrals in the test
-- suite, not just this internal check.
local Complex = require("GABAS-LINAL.modules.complex")

-- Positive abscissas on [-1,1], decreasing, XGK[11] = 0 (the center).
local XGK = {
    0.995657163025808080735527280689003,
    0.973906528517171720077964012084452,
    0.930157491355708226001207180059508,
    0.865063366688984510732096688423493,
    0.780817726586416897063717578345042,
    0.679409568299024406234327365114874,
    0.562757134668604683339000099272694,
    0.433395394129247190799265943165784,
    0.294392862701460198131126603103866,
    0.148874338981631210884826001129720,
    0.000000000000000000000000000000000,
}

-- Kronrod weights, parallel to XGK.
local WGK = {
    0.011694638867371874278064396062192,
    0.032558162307964727478818972459390,
    0.054755896574351996031381300244580,
    0.075039674810919952767043140916190,
    0.093125454583697605535065465083366,
    0.109387158802297641899210590325805,
    0.123491976262065851077958109831074,
    0.134709217311473325928054001771707,
    0.142775938577060080797094273138717,
    0.147739104901338491374841515972068,
    0.149445554002916905664936468389821,
}

-- Gauss weights -- only 5 entries, because the 10-point Gauss rule shares
-- exactly the EVEN-indexed abscissas XGK[2], XGK[4], XGK[6], XGK[8],
-- XGK[10] (the odd-indexed ones, and the center, are Kronrod-only).
local WG = {
    0.066671344308688137593568809893332,
    0.149451349150580593145776339657697,
    0.219086362515982043995534934228163,
    0.269266719309996355091226921569469,
    0.295524224714752870173892994651338,
}

local EPSILON_MACHINE = 2.220446049250313e-16
local UNDERFLOW_THRESHOLD = 2.2250738585072014e-308

local function magnitude(v)
    if Complex.Is_complex(v) then return v:abs() end
    return math.abs(v)
end

-- Claude: evaluates f(x) and classifies the result -- NaN, +-infinity, or
-- a non-numeric return are all reported as a node-evaluation failure
-- (never silently propagated), matching the specification's insistence
-- that "a generic 'integration failed' message is inadequate" (section
-- 19) -- the caller gets told exactly which node and why.
local function safe_eval(f, x)
    local ok, v = pcall(f, x)
    if not ok then
        return nil, tostring(v)
    end
    if Complex.Is_complex(v) then
        if v.re ~= v.re or v.im ~= v.im then
            return nil, "the integrand returned NaN at x=" .. tostring(x)
        end
        if v.re == math.huge or v.re == -math.huge or v.im == math.huge or v.im == -math.huge then
            return nil, "the integrand returned an infinite value at x=" .. tostring(x)
        end
        return v
    end
    if type(v) ~= "number" then
        return nil, "the integrand returned a non-numeric value at x=" .. tostring(x)
    end
    if v ~= v then
        return nil, "the integrand returned NaN at x=" .. tostring(x)
    end
    if v == math.huge or v == -math.huge then
        return nil, "the integrand returned an infinite value at x=" .. tostring(x)
    end
    return v
end

-- Local_Gauss_Kronrod_21(f, left, right): one panel of the specification's
-- own pseudocode (section 9), transcribed directly from QUADPACK's QK21 --
-- same variable names (resg/resk/resabs/resasc/reskh) as the reference,
-- deliberately, so this can be checked line-by-line against the source it
-- is drawn from rather than trusted as an from-scratch reinterpretation.
--
-- Returns: integral_estimate, error_estimate, abs_integral_estimate,
-- abs_deviation_estimate, evaluation_count, error_message (nil on
-- success). half_length is signed (right - left, halved) rather than
-- forced positive, matching the specification's own note that this
-- "reduces the risk of overflow when l and r are large and have the same
-- sign" -- the caller is expected to pass left < right, as
-- Integral_definite's own orientation handling already guarantees.
local function Local_Gauss_Kronrod_21(f, left, right)
    local half_length = (right - left) / 2
    local center = left + half_length

    local fc, fc_err = safe_eval(f, center)
    if fc == nil then
        return nil, nil, nil, nil, 1, fc_err
    end

    local resg = 0
    local resk = WGK[11] * fc
    local resabs = magnitude(resk)
    local evaluations = 1

    local fv1, fv2 = {}, {}

    -- Even-indexed nodes: shared between the Gauss and Kronrod rules.
    for j = 1, 5 do
        local jtw = 2 * j
        local absc = half_length * XGK[jtw]
        local f1, e1 = safe_eval(f, center - absc)
        if f1 == nil then return nil, nil, nil, nil, evaluations, e1 end
        local f2, e2 = safe_eval(f, center + absc)
        if f2 == nil then return nil, nil, nil, nil, evaluations, e2 end
        evaluations = evaluations + 2
        fv1[jtw], fv2[jtw] = f1, f2
        local fsum = f1 + f2
        resg = resg + WG[j] * fsum
        resk = resk + WGK[jtw] * fsum
        resabs = resabs + WGK[jtw] * (magnitude(f1) + magnitude(f2))
    end

    -- Odd-indexed nodes: Kronrod-only, no Gauss contribution.
    for j = 1, 5 do
        local jtwm1 = 2 * j - 1
        local absc = half_length * XGK[jtwm1]
        local f1, e1 = safe_eval(f, center - absc)
        if f1 == nil then return nil, nil, nil, nil, evaluations, e1 end
        local f2, e2 = safe_eval(f, center + absc)
        if f2 == nil then return nil, nil, nil, nil, evaluations, e2 end
        evaluations = evaluations + 2
        fv1[jtwm1], fv2[jtwm1] = f1, f2
        local fsum = f1 + f2
        resk = resk + WGK[jtwm1] * fsum
        resabs = resabs + WGK[jtwm1] * (magnitude(f1) + magnitude(f2))
    end

    local reskh = resk * 0.5
    local resasc = WGK[11] * magnitude(fc - reskh)
    for j = 1, 10 do
        resasc = resasc + WGK[j] * (magnitude(fv1[j] - reskh) + magnitude(fv2[j] - reskh))
    end

    resg = resg * half_length
    resk = resk * half_length
    resabs = resabs * math.abs(half_length)
    resasc = resasc * math.abs(half_length)

    -- Claude: the rescaled-error formula, transcribed from QK21's own
    -- abserr calculation -- NOT just |resk - resg|. A raw rule difference
    -- can be a poor error estimate on its own; this rescaling (and the
    -- roundoff floor) is exactly what the specification's section 6.4
    -- refers to as "considerably safer than blindly using |K-G| as the
    -- final local error."
    local raw_error = magnitude(resk - resg)
    local error = raw_error
    if resasc ~= 0 and raw_error ~= 0 then
        error = resasc * math.min(1, (200 * raw_error / resasc) ^ 1.5)
    end
    if resabs > UNDERFLOW_THRESHOLD / (50 * EPSILON_MACHINE) then
        error = math.max(EPSILON_MACHINE * 50 * resabs, error)
    end

    return resk, error, resabs, resasc, evaluations, nil
end

return {
    Local_Gauss_Kronrod_21 = Local_Gauss_Kronrod_21,
    magnitude = magnitude,
}
