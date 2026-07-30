-- [Vector]
-- Vectors are plain flat tables {v1, v2, ...}, unlike matrices, which are
-- tables of row-tables.
local Core = require("GABAS-LINAL.modules.core")

local function VTrace(matriz)
    local rows, cols = #matriz, #matriz[1]
    assert(rows == cols, "VTrace: matrix must be square.")
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
    assert(#v1 == #v2, "Vdot: vectors must have the same length.")
    local sum = 0
    for i = 1, #v1 do
        sum = sum + v1[i] * v2[i]
    end
    return sum
end

local function VNorm(v, p)
    assert(type(v) == "table" and #v > 0, "VNorm: v must be a non-empty vector.")
    p = p or 2
    assert(p ~= 0, "VNorm: p must not be 0 (there is no p=0 'norm' expressible via sum(|vi|^p)^(1/p)).")
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
    assert(norm > Core.EPSILON, "VNormalize: cannot normalize a zero vector.")
    local C = {}
    for i = 1, #v do
        C[i] = v[i] / norm
    end
    return C
end

return {
    VTrace = VTrace,
    Vdot = Vdot,
    VNorm = VNorm,
    VNormalize = VNormalize,
}
