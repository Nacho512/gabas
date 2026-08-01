-- [Vector]
-- Vectors are plain flat tables {v1, v2, ...}, unlike matrices, which are
-- tables of row-tables.
local Core = require("GABAS-LINAL.modules.core")

local function VTrace(matriz)
    local rows = Core.assert_matrix(matriz, "VTrace: matriz", {square = true})
    local sum = 0
    for i = 1, rows do
        sum = sum + matriz[i][i]
    end
    return sum
end

-- This is the bilinear (non-conjugated) dot product; for the Hermitian
-- inner product on complex vectors, conjugate one side first with
-- Conjugate() (from Complex).
local function Vdot(v1, v2)
    local n = Core.assert_vector(v1, "Vdot: v1")
    Core.assert_vector(v2, "Vdot: v2", {length = n})
    local sum = 0
    for i = 1, #v1 do
        sum = sum + v1[i] * v2[i]
    end
    return sum
end

local function VNorm(v, p)
    assert(type(v) == "table" and #v > 0, "VNorm: v must be a non-empty vector.")
    p = p or 2
    Core.assert_finite_number(p, "VNorm: p")
    assert(p >= 1, "VNorm: p must be at least 1.")
    Core.assert_vector(v, "VNorm: v")
    local sum = 0
    for i = 1, #v do
        sum = sum + Core.cabs(v[i]) ^ p
    end
    return sum ^ (1 / p)
end

-- Returns a new vector scaled to unit length (VNorm(v, p) == 1), leaving `v`
-- itself untouched. Uses the same p as VNorm (default 2, Euclidean).
local function VNormalize(v, p)
    local norm = VNorm(v, p)
    assert(norm > 0, "VNormalize: cannot normalize a zero vector.")
    local C = {}
    for i = 1, #v do
        C[i] = v[i] / norm
    end
    return C
end

-- Cross product: only defined for two 3-dimensional vectors -- the unique
-- vector orthogonal to both v1 and v2, with magnitude ||v1||*||v2||*sin(theta)
-- and orientation given by the right-hand rule. The formula is pure
-- +/-/* on entries, so it works unchanged for Complex-valued vectors too,
-- same as every other function in this module.
local function Cross_product(v1, v2)
    assert(type(v1) == "table" and #v1 == 3, "Cross_product: v1 must be a 3-dimensional vector.")
    assert(type(v2) == "table" and #v2 == 3, "Cross_product: v2 must be a 3-dimensional vector.")
    for i = 1, 3 do
        Core.assert_finite_scalar(v1[i], "Cross_product: v1[" .. i .. "]")
        Core.assert_finite_scalar(v2[i], "Cross_product: v2[" .. i .. "]")
    end
    return {
        v1[2] * v2[3] - v1[3] * v2[2],
        v1[3] * v2[1] - v1[1] * v2[3],
        v1[1] * v2[2] - v1[2] * v2[1],
    }
end

return {
    V_trace = VTrace,
    V_dot = Vdot,
    V_norm = VNorm,
    V_normalize = VNormalize,
    Cross_product = Cross_product,
}
