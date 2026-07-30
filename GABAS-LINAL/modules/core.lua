-- [Core]
-- Internal infrastructure shared by 2+ of the domain submodules (Matrix,
-- Tensor, Vector, Determinant, Systems, Eigen, FFT). Nothing in this file
-- is meant to be part of GABAS_LINAL.lua's public API -- every function
-- here is plumbing, not something an end user calls directly, which is why
-- GABAS_LINAL.lua's aggregator loop deliberately excludes "core" from the
-- list of submodules it flattens into the public table.
local Complex = require("GABAS-LINAL.modules.complex")

local EPSILON = 1e-9

-- Magnitude that works on both plain numbers and Complex values -- used
-- everywhere a pivoting/convergence/near-zero check needs "how big is
-- this", regardless of whether the entries are real or complex.
local function cabs(x)
    if Complex.Is_complex(x) then
        return x:abs()
    end
    return math.abs(x)
end

-- Used by MDet/Big_MDet (Determinant), Rank (Systems), and
-- Eigenvalues/Eigenvectors (Eigen) alike, so callers never mutate the
-- matrix they passed in.
local function copy_matrix(matriz)
    local C = {}
    for i = 1, #matriz do
        local row = {}
        for j = 1, #matriz[i] do
            row[j] = matriz[i][j]
        end
        C[i] = row
    end
    return C
end

-- Searches rows [k, last_row] of column `col` for the entry with the
-- largest magnitude (classic partial pivoting: picking the largest
-- magnitude, not just the first non-zero one, keeps the elimination
-- numerically stable). Returns the row index, or nil if every candidate is
-- ~0 (singular column). Shared by Determinant and Systems.
local function find_pivot_row(mat, k, last_row, col)
    local best_row, best_val = nil, EPSILON
    for i = k, last_row do
        local v = cabs(mat[i][col])
        if v > best_val then
            best_row, best_val = i, v
        end
    end
    return best_row
end

-- Shared by Big_MDet (Determinant) and FFT -- a genuine cross-category
-- dependency, which is exactly why next_pow2 lives here instead of inside
-- either of those two submodules.
local function next_pow2(x)
    local m = 1
    while m < x do m = m * 2 end
    return m
end

return {
    EPSILON = EPSILON,
    cabs = cabs,
    copy_matrix = copy_matrix,
    find_pivot_row = find_pivot_row,
    next_pow2 = next_pow2,
}
