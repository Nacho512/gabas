-- [Eigen]
local Core = require("GABAS-LINAL.modules.core")
local Complex = require("GABAS-LINAL.modules.complex")
local Matrix = require("GABAS-LINAL.modules.matrix")
local Vector = require("GABAS-LINAL.modules.vector")

-- Modified Gram-Schmidt QR decomposition of a REAL square matrix (Q
-- orthonormal, R upper triangular). Internal helper for Eigenvalues/
-- Eigenvectors. Deliberately real-only: correct QR for complex matrices
-- needs a conjugated ("Hermitian") Gram-Schmidt, which isn't implemented
-- here -- Eigenvalues/Eigenvectors reject Complex-entry input accordingly.
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
        assert(norm > Core.EPSILON, "qr_decompose: matrix is numerically singular.")
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

local function assert_real_matrix(matriz, fn_name)
    for i = 1, #matriz do
        for j = 1, #matriz[i] do
            assert(not Complex.Is_complex(matriz[i][j]),
                fn_name .. ": Complex-entry input matrices are not supported (the QR " ..
                "iteration only runs in real arithmetic); Complex eigenvalues can still " ..
                "occur in the *output* for a real matrix with a complex-conjugate pair.")
        end
    end
end

-- Unshifted QR algorithm: iterating A <- R*Q (from A's QR decomposition)
-- converges A towards a quasi-triangular (real Schur) form. Any diagonal
-- entry whose subdiagonal neighbor vanished is a converged real eigenvalue.
-- A 2x2 block that refuses to converge encodes a complex-conjugate
-- eigenvalue pair; its two roots are extracted analytically via the
-- quadratic formula on that block's trace/determinant (the same technique
-- LAPACK-style real-Schur eigensolvers use), which is exact regardless of
-- how far the QR iteration itself has converged. This is what actually
-- resolves the old "won't converge for complex eigenvalues" limitation,
-- without needing a complex-arithmetic QR implementation.
local function Eigenvalues(matriz, max_iter, tol)
    local n = #matriz
    assert(n == #matriz[1], "Eigenvalues: matrix must be square.")
    assert_real_matrix(matriz, "Eigenvalues")
    max_iter = max_iter or 500
    -- Silently accepting max_iter <= 0 would just skip the whole iteration
    -- and return the unconverged input read straight off the diagonal --
    -- wrong, but without any error to say so.
    assert(type(max_iter) == "number" and max_iter >= 1 and max_iter == math.floor(max_iter),
        "Eigenvalues: max_iter must be a positive integer.")
    tol = tol or 1e-10

    local A = Core.copy_matrix(matriz)

    for _ = 1, max_iter do
        local Q, R = qr_decompose(A)
        A = Matrix.Mat_mul(R, Q)

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
    local i = 1
    while i <= n do
        if i < n and math.abs(A[i + 1][i]) > tol then
            local a, b = A[i][i], A[i][i + 1]
            local c, d = A[i + 1][i], A[i + 1][i + 1]
            local tr = a + d
            local det = a * d - b * c
            local disc = tr * tr - 4 * det
            if disc >= 0 then
                local sq = math.sqrt(disc)
                eigenvalues[i] = (tr + sq) / 2
                eigenvalues[i + 1] = (tr - sq) / 2
            else
                local sq = math.sqrt(-disc)
                eigenvalues[i] = Complex.Complex(tr / 2, sq / 2)
                eigenvalues[i + 1] = Complex.Complex(tr / 2, -sq / 2)
            end
            i = i + 2
        else
            eigenvalues[i] = A[i][i]
            i = i + 1
        end
    end
    return eigenvalues
end

-- Same iteration as Eigenvalues, additionally accumulating the orthogonal
-- factors Q_total = Q_1 * Q_2 * ... ; its columns converge to the
-- eigenvectors. Rigorously guaranteed for symmetric matrices; treat the
-- vectors as approximate for general (non-symmetric) matrices. For a matrix
-- with a complex-conjugate eigenvalue pair, the corresponding two columns
-- of Q_total only span the associated real 2-D invariant subspace, not the
-- individual complex eigenvectors -- use Eigenvalues() for the (possibly
-- complex) eigenvalues themselves; extracting the true complex eigenvectors
-- for those pairs isn't implemented yet. Returns eigenvalues, and a matrix
-- whose column k is the eigenvector for eigenvalues[k].
local function Eigenvectors(matriz, max_iter, tol)
    local n = #matriz
    assert(n == #matriz[1], "Eigenvectors: matrix must be square.")
    assert_real_matrix(matriz, "Eigenvectors")
    max_iter = max_iter or 500
    assert(type(max_iter) == "number" and max_iter >= 1 and max_iter == math.floor(max_iter),
        "Eigenvectors: max_iter must be a positive integer.")
    tol = tol or 1e-10

    local A = Core.copy_matrix(matriz)
    local Q_total = Matrix.Eye(n)

    for _ = 1, max_iter do
        local Q, R = qr_decompose(A)
        A = Matrix.Mat_mul(R, Q)
        Q_total = Matrix.Mat_mul(Q_total, Q)

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

-- Modified Gram-Schmidt: takes a matrix whose ROWS are a (possibly linearly
-- dependent, possibly Complex-valued) set of input vectors, and returns a
-- matrix whose rows are an orthonormal basis for their span -- using the
-- Hermitian inner product (Conjugate(qi) . v) for the projection
-- coefficients, so this works correctly for Complex-entry input too, not
-- just real. A vector that turns out to be a linear combination of the ones
-- already processed (residual norm below `tol` after projecting out every
-- prior basis vector) is skipped rather than raising an error -- the result
-- can have fewer rows than the input; a warning is printed (and the count
-- returned) when that happens, so it's visible rather than silently losing
-- a vector.
local function GRAM_SCH(vectores, tol)
    assert(type(vectores) == "table" and #vectores > 0, "GRAM_SCH: need at least one vector.")
    local k = #vectores
    assert(type(vectores[1]) == "table" and #vectores[1] > 0, "GRAM_SCH: vectors must be non-empty.")
    local dim = #vectores[1]
    for i = 2, k do
        assert(type(vectores[i]) == "table" and #vectores[i] == dim,
            "GRAM_SCH: all vectors must have the same length.")
    end
    tol = tol or 1e-10

    local basis = {}
    for i = 1, k do
        local v = {}
        for d = 1, dim do
            v[d] = vectores[i][d]
        end
        for _, q in ipairs(basis) do
            local coeff = Vector.Vdot(Complex.Conjugate(q), v)
            for d = 1, dim do
                v[d] = v[d] - coeff * q[d]
            end
        end
        if Vector.VNorm(v) > tol then
            basis[#basis + 1] = Vector.VNormalize(v)
        end
    end

    local skipped = k - #basis
    assert(#basis > 0, "GRAM_SCH: all input vectors were linearly dependent (zero-dimensional span).")
    if skipped > 0 then
        io.stderr:write(string.format(
            "GRAM_SCH: %d of %d input vector(s) were linearly dependent on the others and were dropped.\n",
            skipped, k))
    end

    return basis, skipped
end

return {
    Eigenvalues = Eigenvalues,
    Eigenvectors = Eigenvectors,
    GRAM_SCH = GRAM_SCH,
}
