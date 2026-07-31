-- [Determinant]
local Core = require("GABAS-LINAL.modules.core")
local Matrix = require("GABAS-LINAL.modules.matrix")

local function MDet(matrix)
    local rows = Core.assert_matrix(matrix, "MDet: matrix", {square = true})

    -- LU decomposition (Gaussian elimination) with partial pivoting.
    local lu = Core.copy_matrix(matrix)
    local sign = 1

    for k = 1, rows - 1 do
        local pivot_row = Core.find_pivot_row(lu, k, rows, k)
        if not pivot_row then
            -- Entire column below (and at) the pivot is zero: matrix is singular.
            return 0
        end
        if pivot_row ~= k then
            lu[k], lu[pivot_row] = lu[pivot_row], lu[k]
            sign = -sign
        end

        for i = k + 1, rows do
            local factor = lu[i][k] / lu[k][k]
            for j = k + 1, rows do
                lu[i][j] = lu[i][j] - factor * lu[k][j]
            end
            lu[i][k] = factor
        end
    end

    if Core.cabs(lu[rows][rows]) == 0 then
        return 0
    end

    local det = sign
    for i = 1, rows do
        det = det * lu[i][i]
    end

    return det
end

-- ===========================================================================
-- Big_MDet: determinant of huge matrices via blocked LU decomposition with
-- partial pivoting, where the dominant O(n^3) step -- updating the trailing
-- submatrix after each panel -- is computed with Strassen's algorithm
-- (~O(n^2.807)) instead of naive triple-loop multiplication. This is the
-- standard "GEPP + Strassen" hybrid design real numerical libraries use: full
-- panels are still factored with ordinary partial pivoting (so it's exactly
-- as numerically stable as MDet -- there is no bad-pivot fragility here),
-- Strassen only accelerates the large rank-`block_size` update, which is
-- where virtually all the arithmetic actually happens for big n.
--
-- Below `block_size`, and for the recursion's own internal base case, this
-- just falls back to MDet's plain pivoted elimination / Mat_mul -- Strassen's
-- overhead (recursive submatrix allocation, 7 sub-multiplies + ~18 add/sub
-- passes per level) only pays off once n is large enough to amortize it.
-- Caveat worth being upfront about: in pure interpreted Lua, table-of-tables
-- allocation at every recursion level carries real constant-factor cost, so
-- the practical crossover where this beats MDet outright is workload- and
-- machine-dependent -- likely somewhere in the hundreds-to-low-thousands of
-- rows, not the dozens where Strassen starts winning in a compiled, cache-
-- friendly implementation. Tune `block_size` (default 64) if needed.
-- ===========================================================================

-- Square Strassen multiply; n must already be a power of 2.
local function strassen_mul(A, B, n, threshold)
    if n <= threshold then
        return Matrix.Mat_mul(A, B)
    end
    local h = n // 2
    local A11, A12, A21, A22 = {}, {}, {}, {}
    local B11, B12, B21, B22 = {}, {}, {}, {}
    for i = 1, h do
        A11[i], A12[i], A21[i], A22[i] = {}, {}, {}, {}
        B11[i], B12[i], B21[i], B22[i] = {}, {}, {}, {}
        for j = 1, h do
            A11[i][j] = A[i][j]
            A12[i][j] = A[i][j + h]
            A21[i][j] = A[i + h][j]
            A22[i][j] = A[i + h][j + h]
            B11[i][j] = B[i][j]
            B12[i][j] = B[i][j + h]
            B21[i][j] = B[i + h][j]
            B22[i][j] = B[i + h][j + h]
        end
    end

    local M1 = strassen_mul(Matrix.Mat_sum(A11, A22), Matrix.Mat_sum(B11, B22), h, threshold)
    local M2 = strassen_mul(Matrix.Mat_sum(A21, A22), B11, h, threshold)
    local M3 = strassen_mul(A11, Matrix.Mat_sub(B12, B22), h, threshold)
    local M4 = strassen_mul(A22, Matrix.Mat_sub(B21, B11), h, threshold)
    local M5 = strassen_mul(Matrix.Mat_sum(A11, A12), B22, h, threshold)
    local M6 = strassen_mul(Matrix.Mat_sub(A21, A11), Matrix.Mat_sum(B11, B12), h, threshold)
    local M7 = strassen_mul(Matrix.Mat_sub(A12, A22), Matrix.Mat_sum(B21, B22), h, threshold)

    local C11 = Matrix.Mat_sum(Matrix.Mat_sub(Matrix.Mat_sum(M1, M4), M5), M7)
    local C12 = Matrix.Mat_sum(M3, M5)
    local C21 = Matrix.Mat_sum(M2, M4)
    local C22 = Matrix.Mat_sum(Matrix.Mat_sub(Matrix.Mat_sum(M1, M3), M2), M6)

    local C = {}
    for i = 1, h do
        local row_top, row_bot = {}, {}
        for j = 1, h do
            row_top[j], row_top[j + h] = C11[i][j], C12[i][j]
            row_bot[j], row_bot[j + h] = C21[i][j], C22[i][j]
        end
        C[i], C[i + h] = row_top, row_bot
    end
    return C
end

-- Rectangular multiply (m x k) * (k x n) via square Strassen, using
-- zero-padding up to the next power of 2 and extracting the real m x n
-- corner of the result -- the zero-padded rows/columns only ever contribute
-- zero to the extracted region, so this is exact, not approximate.
local function strassen_mul_rect(A, m, k, B, k2, n, threshold)
    assert(k == k2, "strassen_mul_rect: inner dimensions must match.")
    local size = Core.next_pow2(math.max(m, k, n))
    local Ap, Bp = {}, {}
    for i = 1, size do
        local row = {}
        for j = 1, size do
            row[j] = (i <= m and j <= k) and A[i][j] or 0
        end
        Ap[i] = row
    end
    for i = 1, size do
        local row = {}
        for j = 1, size do
            row[j] = (i <= k and j <= n) and B[i][j] or 0
        end
        Bp[i] = row
    end
    local Cp = strassen_mul(Ap, Bp, size, threshold)
    local C = {}
    for i = 1, m do
        local row = {}
        for j = 1, n do
            row[j] = Cp[i][j]
        end
        C[i] = row
    end
    return C
end

local function Big_MDet(matrix, block_size)
    local n = Core.assert_matrix(matrix, "Big_MDet: matrix", {square = true})
    block_size = block_size or 64
    -- A block_size <= 0 would make panel_end < col every iteration, so `col`
    -- never advances -- an infinite loop, not just a bad result.
    Core.assert_positive_integer(block_size, "Big_MDet: block_size")

    if n <= block_size then
        return MDet(matrix)
    end

    local lu = Core.copy_matrix(matrix)
    local sign = 1
    local col = 1

    while col <= n do
        local panel_end = math.min(col + block_size - 1, n)

        -- Panel factorization: ordinary partial-pivoting elimination,
        -- touching only this panel's own columns (col..panel_end) for every
        -- row, panel or trailing alike -- deliberately NOT the full row
        -- width. This has to stay column-restricted here because a pivot
        -- swap on a *later* k within this same panel can still move a
        -- "trailing" row into panel territory; updating beyond-panel columns
        -- eagerly, keyed off a row's position at the time, silently uses the
        -- wrong (pre-swap) data for whichever row a later swap relocates.
        -- Row swaps themselves are still full-row (so A12 ends up correctly
        -- permuted), just not the column *updates*.
        for k = col, panel_end do
            local pivot_row = Core.find_pivot_row(lu, k, n, k)
            if not pivot_row then
                return 0
            end
            if pivot_row ~= k then
                lu[k], lu[pivot_row] = lu[pivot_row], lu[k]
                sign = -sign
            end
            for i = k + 1, n do
                local factor = lu[i][k] / lu[k][k]
                for j = k + 1, panel_end do
                    lu[i][j] = lu[i][j] - factor * lu[k][j]
                end
                lu[i][k] = factor
            end
        end

        -- U12 = L11^-1 * A12, computed now that every row's final position
        -- within the panel is settled (no more swaps will touch rows
        -- col..panel_end after this point). A12 -- panel rows, columns
        -- beyond panel_end -- was left untouched above, so it still holds
        -- the original (but already correctly row-permuted) values; replay
        -- the panel's already-stored multipliers against just those columns,
        -- restricted to the `width`-many panel rows -- cheap regardless of n.
        for k = col, panel_end - 1 do
            for i = k + 1, panel_end do
                local factor = lu[i][k]
                for j = panel_end + 1, n do
                    lu[i][j] = lu[i][j] - factor * lu[k][j]
                end
            end
        end

        -- Trailing-submatrix update: instead of continuing the panel's
        -- row-by-row elimination out to column n (the O(n^3)-dominant part
        -- of plain Gaussian elimination), do it as one big
        -- Strassen-accelerated matrix multiply.
        if panel_end < n then
            local rest = n - panel_end
            local width = panel_end - col + 1
            local L21 = {}
            for i = 1, rest do
                local row = {}
                for j = 1, width do
                    row[j] = lu[panel_end + i][col + j - 1]
                end
                L21[i] = row
            end
            local U12 = {}
            for i = 1, width do
                local row = {}
                for j = 1, rest do
                    row[j] = lu[col + i - 1][panel_end + j]
                end
                U12[i] = row
            end
            local update = strassen_mul_rect(L21, rest, width, U12, width, rest, block_size)
            for i = 1, rest do
                local trailing_row = lu[panel_end + i]
                local update_row = update[i]
                for j = 1, rest do
                    trailing_row[panel_end + j] = trailing_row[panel_end + j] - update_row[j]
                end
            end
        end

        col = panel_end + 1
    end

    if Core.cabs(lu[n][n]) == 0 then
        return 0
    end

    local det = sign
    for i = 1, n do
        det = det * lu[i][i]
    end

    return det
end

return {
    MDet = MDet,
    M_Det = MDet,
    Big_MDet = Big_MDet,
    Big_M_Det = Big_MDet,
}
