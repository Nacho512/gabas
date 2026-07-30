local EPSILON = 1e-9

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

-- Searches rows [k, last_row] of column `col` for the entry with the largest
-- absolute value (classic partial pivoting: picking the largest magnitude,
-- not just the first non-zero one, keeps the elimination numerically stable).
-- Returns the row index, or nil if every candidate is ~0 (singular column).
local function find_pivot_row(mat, k, last_row, col)
    local best_row, best_val = nil, EPSILON
    for i = k, last_row do
        local v = math.abs(mat[i][col])
        if v > best_val then
            best_row, best_val = i, v
        end
    end
    return best_row
end

local function Show_matrix(matriz)
    local m, n = #matriz, #matriz[1]
    for i = 1, m do
        local row = matriz[i]
        for j = 1, n do
            io.write(row[j], "\t")
        end
        print()
    end
end

local function Matrix(data, num_filas)
    local limite = #data
    assert((limite % num_filas) == 0, "Mismatch between number of elements and number of rows.")
    local m = {}
    local prop = limite / num_filas
    local k = 1
    for i = 1, num_filas do
        m[i] = {}
        for j = 1, prop do
            m[i][j] = data[k]
            k = k + 1
        end
    end
    return m
end

local function Sequence(aleph, tat)
    local m = {}
    local k = 1
    for j = aleph, tat do
        m[k] = j
        k = k + 1
    end
    return m
end

local function Show_table(tabla)
    local limite = #tabla
    for j = 1, limite do
        print(tabla[j])
    end
end

local function Mat_mul(M1, M2)
    local m, p = #M1, #M1[1]
    local n = #M2[1]
    assert(#M2 == p, "Mat_mul: number of columns in M1 must match number of rows in M2.")
    local C = {}
    for i = 1, m do
        local row1 = M1[i]
        local row_c = {}
        for j = 1, n do
            local sum = 0
            for k = 1, p do
                sum = sum + row1[k] * M2[k][j]
            end
            row_c[j] = sum
        end
        C[i] = row_c
    end
    return C
end

local function Mat_sum(M1, M2)
    local m, n = #M1, #M1[1]
    assert(#M2 == m and #M2[1] == n, "Mat_sum: matrices must have matching dimensions.")
    local C = {}
    for i = 1, m do
        local row1, row2 = M1[i], M2[i]
        local row_c = {}
        for j = 1, n do
            row_c[j] = row1[j] + row2[j]
        end
        C[i] = row_c
    end
    return C
end

local function Mat_sub(M1, M2)
    local m, n = #M1, #M1[1]
    assert(#M2 == m and #M2[1] == n, "Mat_sub: matrices must have matching dimensions.")
    local C = {}
    for i = 1, m do
        local row1, row2 = M1[i], M2[i]
        local row_c = {}
        for j = 1, n do
            row_c[j] = row1[j] - row2[j]
        end
        C[i] = row_c
    end
    return C
end

local function Scalar_mul(matriz, escalar)
    local m, n = #matriz, #matriz[1]
    local C = {}
    for i = 1, m do
        local row = matriz[i]
        local row_c = {}
        for j = 1, n do
            row_c[j] = row[j] * escalar
        end
        C[i] = row_c
    end
    return C
end

local function Eye(dim)
    local C = {}
    for i = 1, dim do
        C[i] = {}
        for j = 1, dim do
            if i == j then
                C[i][j] = 1
            else
                C[i][j] = 0
            end
        end
    end
    return C
end

local function Zeroes(n_rows, n_cols)
    local C = {}
    for i = 1, n_rows do
        C[i] = {}
        for j = 1, n_cols do
            C[i][j] = 0
        end
    end
    return C
end

local function Random_mat(m, n, valor)
    local C = {}
    for i = 1, m do
        C[i] = {}
        for j = 1, n do
            C[i][j] = math.random(1, valor)
        end
    end
    return C
end

local function T(matriz)
    local m, n = #matriz, #matriz[1]
    local C = {}
    for i = 1, n do
        local row_c = {}
        for j = 1, m do
            row_c[j] = matriz[j][i]
        end
        C[i] = row_c
    end
    return C
end

local function Trace(matriz)
    local rows, cols = #matriz, #matriz[1]
    assert(rows == cols, "Trace: matrix must be square.")
    local sum = 0
    for i = 1, rows do
        sum = sum + matriz[i][i]
    end
    return sum
end

-- Vectors are plain flat tables {v1, v2, ...}, unlike matrices, which are
-- tables of row-tables.
local function Dot(v1, v2)
    assert(#v1 == #v2, "Dot: vectors must have the same length.")
    local sum = 0
    for i = 1, #v1 do
        sum = sum + v1[i] * v2[i]
    end
    return sum
end

local function Norm(v, p)
    p = p or 2
    local sum = 0
    for i = 1, #v do
        sum = sum + math.abs(v[i]) ^ p
    end
    return sum ^ (1 / p)
end

local function Determinant(matrix)
    local rows, cols = #matrix, #matrix[1]
    assert(rows == cols, "Determinant: matrix must be square.")

    -- LU decomposition (Gaussian elimination) with partial pivoting.
    local lu = copy_matrix(matrix)
    local sign = 1

    for k = 1, rows - 1 do
        local pivot_row = find_pivot_row(lu, k, rows, k)
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

    if math.abs(lu[rows][rows]) <= EPSILON then
        return 0
    end

    local det = sign
    for i = 1, rows do
        det = det * lu[i][i]
    end

    return det
end

-- Row-reduces a copy of the matrix (Gaussian elimination with partial
-- pivoting) and counts the pivots found. Works for non-square matrices too.
local function Rank(matriz)
    local mat = copy_matrix(matriz)
    local rows, cols = #mat, #mat[1]
    local rank = 0

    for col = 1, cols do
        if rank >= rows then break end
        local pivot_row = find_pivot_row(mat, rank + 1, rows, col)
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
local function Inverse(matriz)
    local n = #matriz
    assert(n == #matriz[1], "Inverse: matrix must be square.")

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
        local pivot_row = find_pivot_row(aug, k, n, k)
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
                if factor ~= 0 then
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
    local n = #A
    assert(n == #A[1], "Solve: coefficient matrix A must be square.")
    assert(#b == n, "Solve: vector b must have the same length as the number of rows in A.")

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
        local pivot_row = find_pivot_row(aug, k, n, k)
        assert(pivot_row, "Solve: matrix A is singular; the system has no unique solution.")
        if pivot_row ~= k then
            aug[k], aug[pivot_row] = aug[pivot_row], aug[k]
        end

        for i = k + 1, n do
            local factor = aug[i][k] / aug[k][k]
            if factor ~= 0 then
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

-- Modified Gram-Schmidt QR decomposition of a square matrix (Q orthonormal,
-- R upper triangular). Internal helper for Eigenvalues/Eigenvectors.
local function qr_decompose(mat)
    local n = #mat

    local q_cols = {}
    for j = 1, n do
        local v = {}
        for i = 1, n do
            v[i] = mat[i][j]
        end
        q_cols[j] = v
    end

    local R = {}
    for i = 1, n do
        R[i] = {}
        for j = 1, n do
            R[i][j] = 0
        end
    end

    for j = 1, n do
        for i = 1, j - 1 do
            local dot = 0
            for r = 1, n do
                dot = dot + q_cols[i][r] * q_cols[j][r]
            end
            R[i][j] = dot
            for r = 1, n do
                q_cols[j][r] = q_cols[j][r] - dot * q_cols[i][r]
            end
        end
        local norm = 0
        for r = 1, n do
            norm = norm + q_cols[j][r] ^ 2
        end
        norm = math.sqrt(norm)
        assert(norm > EPSILON, "qr_decompose: matrix is numerically singular.")
        R[j][j] = norm
        for r = 1, n do
            q_cols[j][r] = q_cols[j][r] / norm
        end
    end

    local Q = {}
    for i = 1, n do
        local row = {}
        for j = 1, n do
            row[j] = q_cols[j][i]
        end
        Q[i] = row
    end

    return Q, R
end

-- Unshifted QR algorithm: iterating A <- R*Q (from A's QR decomposition)
-- converges the diagonal of A to its eigenvalues. This only works for real
-- eigenvalues -- there is no complex-number support here, so a matrix with
-- complex-conjugate eigenvalue pairs (rotation-like matrices, for instance)
-- will not converge and the returned "eigenvalues" will be inaccurate. It's
-- reliable for symmetric matrices and any matrix whose eigenvalues are all
-- real and of distinct magnitude.
local function Eigenvalues(matriz, max_iter, tol)
    local n = #matriz
    assert(n == #matriz[1], "Eigenvalues: matrix must be square.")
    max_iter = max_iter or 500
    tol = tol or 1e-10

    local A = copy_matrix(matriz)

    for _ = 1, max_iter do
        local Q, R = qr_decompose(A)
        A = Mat_mul(R, Q)

        local off_diag = 0
        for i = 2, n do
            for j = 1, i - 1 do
                off_diag = off_diag + math.abs(A[i][j])
            end
        end
        if off_diag < tol then
            break
        end
    end

    local eigenvalues = {}
    for i = 1, n do
        eigenvalues[i] = A[i][i]
    end
    return eigenvalues
end

-- Same iteration as Eigenvalues, additionally accumulating the orthogonal
-- factors Q_total = Q_1 * Q_2 * ... ; its columns converge to the
-- eigenvectors. Rigorously guaranteed for symmetric matrices; treat the
-- vectors as approximate for general (non-symmetric) matrices. Returns
-- eigenvalues, and a matrix whose column k is the eigenvector for
-- eigenvalues[k].
local function Eigenvectors(matriz, max_iter, tol)
    local n = #matriz
    assert(n == #matriz[1], "Eigenvectors: matrix must be square.")
    max_iter = max_iter or 500
    tol = tol or 1e-10

    local A = copy_matrix(matriz)
    local Q_total = Eye(n)

    for _ = 1, max_iter do
        local Q, R = qr_decompose(A)
        A = Mat_mul(R, Q)
        Q_total = Mat_mul(Q_total, Q)

        local off_diag = 0
        for i = 2, n do
            for j = 1, i - 1 do
                off_diag = off_diag + math.abs(A[i][j])
            end
        end
        if off_diag < tol then
            break
        end
    end

    local eigenvalues = {}
    for i = 1, n do
        eigenvalues[i] = A[i][i]
    end

    return eigenvalues, Q_total
end

return {
    Show_matrix = Show_matrix, Matrix = Matrix, Sequence = Sequence, Show_table = Show_table,
    Mat_mul = Mat_mul, Mat_sum = Mat_sum, Mat_sub = Mat_sub, Scalar_mul = Scalar_mul,
    Eye = Eye, Zeroes = Zeroes, Random_mat = Random_mat, T = T, Trace = Trace,
    Dot = Dot, Norm = Norm, Determinant = Determinant, Rank = Rank, Inverse = Inverse,
    Solve = Solve, Eigenvalues = Eigenvalues, Eigenvectors = Eigenvectors,
}
