-- [Systems]
local Core = require("GABAS-LINAL.modules.core")

-- Row-reduces a copy of the matrix (Gaussian elimination with partial
-- pivoting) and counts the pivots found. Works for non-square matrices too.
local function Rank(matriz)
    Core.assert_matrix(matriz, "Rank: matriz")
    local mat = Core.copy_matrix(matriz)
    local rows, cols = #mat, #mat[1]
    local rank = 0

    for col = 1, cols do
        if rank >= rows then break end
        local pivot_row = Core.find_pivot_row(mat, rank + 1, rows, col)
        if pivot_row then
            if pivot_row ~= rank + 1 then
                mat[rank + 1], mat[pivot_row] = mat[pivot_row], mat[rank + 1]
            end
            rank = rank + 1
            for i = rank + 1, rows do
                local factor = mat[i][col] / mat[rank][col]
                for j = col, cols do
                    mat[i][j] = mat[i][j] - factor * mat[rank][j]
                end
            end
        end
    end

    return rank
end

-- Gauss-Jordan elimination on [matrix | identity] with partial pivoting.
-- Complex entries work unchanged (division/subtraction below flow through
-- Complex's own operator overloads); pivoting uses Core.cabs, which
-- already handles both real and Complex magnitudes.
local function Inverse(matriz)
    local n = Core.assert_matrix(matriz, "Inverse: matriz", {square = true})

    local aug = {}
    for i = 1, n do
        local row = {}
        for j = 1, n do
            row[j] = matriz[i][j]
        end
        for j = 1, n do
            row[n + j] = (i == j) and 1 or 0
        end
        aug[i] = row
    end

    for k = 1, n do
        local pivot_row = Core.find_pivot_row(aug, k, n, k)
        assert(pivot_row, "Inverse: matrix is singular and cannot be inverted.")
        if pivot_row ~= k then
            aug[k], aug[pivot_row] = aug[pivot_row], aug[k]
        end

        local pivot_row_data = aug[k]
        local pivot = pivot_row_data[k]
        for j = k, 2 * n do
            pivot_row_data[j] = pivot_row_data[j] / pivot
        end

        for i = 1, n do
            if i ~= k then
                local row_i = aug[i]
                local factor = row_i[k]
                if not Core.is_near_zero(factor, Core.cabs(factor), 0, Core.MACHINE_EPSILON) then
                    for j = k, 2 * n do
                        row_i[j] = row_i[j] - factor * pivot_row_data[j]
                    end
                end
            end
        end
    end

    local C = {}
    for i = 1, n do
        local row_c = {}
        for j = 1, n do
            row_c[j] = aug[i][n + j]
        end
        C[i] = row_c
    end
    return C
end

-- Solves Ax = b via Gaussian elimination with partial pivoting and back
-- substitution: O(n^3) and numerically far more stable than expanding
-- Cramer's rule for anything beyond 2x2/3x3 systems, so that's what's used
-- here. `b` is a flat vector, not a matrix.
local function Solve(A, b)
    local n = Core.assert_matrix(A, "Solve: coefficient matrix A", {square = true})
    Core.assert_vector(b, "Solve: vector b", {length = n})

    local aug = {}
    for i = 1, n do
        local row = {}
        for j = 1, n do
            row[j] = A[i][j]
        end
        row[n + 1] = b[i]
        aug[i] = row
    end

    for k = 1, n do
        local pivot_row = Core.find_pivot_row(aug, k, n, k)
        assert(pivot_row, "Solve: matrix A is singular; the system has no unique solution.")
        if pivot_row ~= k then
            aug[k], aug[pivot_row] = aug[pivot_row], aug[k]
        end

        for i = k + 1, n do
            local factor = aug[i][k] / aug[k][k]
            if not Core.is_near_zero(factor, Core.cabs(factor), 0, Core.MACHINE_EPSILON) then
                for j = k, n + 1 do
                    aug[i][j] = aug[i][j] - factor * aug[k][j]
                end
            end
        end
    end

    local x = {}
    for i = n, 1, -1 do
        local sum = aug[i][n + 1]
        for j = i + 1, n do
            sum = sum - aug[i][j] * x[j]
        end
        x[i] = sum / aug[i][i]
    end

    return x
end

return {
    Rank = Rank,
    Inverse = Inverse,
    Solve = Solve,
}
