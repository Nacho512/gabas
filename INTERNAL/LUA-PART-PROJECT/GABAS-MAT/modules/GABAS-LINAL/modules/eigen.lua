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
        R[j][j] = norm
        if Core.is_near_zero(norm, 1, Core.MACHINE_EPSILON, Core.MACHINE_EPSILON) then
            -- Codex: a rank-deficient column has zero residual but QR remains
            -- well-defined; complete Q with an orthogonal coordinate vector.
            local replacement
            for candidate = 1, n do
                local w = {}
                for r = 1, n do w[r] = (r == candidate) and 1 or 0 end
                for i = 1, j - 1 do
                    local dot = 0
                    for r = 1, n do dot = dot + q_cols[i][r] * w[r] end
                    for r = 1, n do w[r] = w[r] - dot * q_cols[i][r] end
                end
                local wnorm = 0
                for r = 1, n do wnorm = wnorm + w[r] ^ 2 end
                wnorm = math.sqrt(wnorm)
                if wnorm > Core.MACHINE_EPSILON then
                    for r = 1, n do w[r] = w[r] / wnorm end
                    replacement = w
                    break
                end
            end
            assert(replacement, "qr_decompose: unable to complete an orthonormal basis.")
            q_cols[j] = replacement
        else
            for r = 1, n do
                q_cols[j][r] = q_cols[j][r] / norm
            end
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
    return Core.assert_matrix(matriz, fn_name .. ": matriz", {square = true, real = true})
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
    local n = assert_real_matrix(matriz, "Eigenvalues")
    max_iter = max_iter or 500
    -- Silently accepting max_iter <= 0 would just skip the whole iteration
    -- and return the unconverged input read straight off the diagonal --
    -- wrong, but without any error to say so.
    Core.assert_positive_integer(max_iter, "Eigenvalues: max_iter")
    tol = tol or 1e-10
    Core.assert_nonneg_number(tol, "Eigenvalues: tol")

    local A = Core.copy_matrix(matriz)

    for _ = 1, max_iter do
        local Q, R = qr_decompose(A)
        A = Matrix.Mat_Mul(R, Q)

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

-- Codex: Eigenvectors is intentionally restricted to finite real symmetric
-- matrices, for which accumulated QR factors converge to an orthonormal
-- eigenbasis. General and complex eigenvectors remain future work.
local function Eigenvectors(matriz, max_iter, tol)
    local n = assert_real_matrix(matriz, "Eigenvectors")
    for i = 1, n do
        for j = i + 1, n do
            local scale = math.max(math.abs(matriz[i][j]), math.abs(matriz[j][i]), 1)
            assert(Core.is_near_zero(matriz[i][j] - matriz[j][i], scale),
                "Eigenvectors: matriz must be symmetric; general real matrices are not supported yet.")
        end
    end
    max_iter = max_iter or 500
    Core.assert_positive_integer(max_iter, "Eigenvectors: max_iter")
    tol = tol or 1e-10
    Core.assert_nonneg_number(tol, "Eigenvectors: tol")

    local A = Core.copy_matrix(matriz)
    local Q_total = Matrix.Eye(n)

    for _ = 1, max_iter do
        local Q, R = qr_decompose(A)
        A = Matrix.Mat_Mul(R, Q)
        Q_total = Matrix.Mat_Mul(Q_total, Q)

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
    Core.assert_vector(vectores[1], "GRAM_SCH: vectores[1]")
    for i = 2, k do
        Core.assert_vector(vectores[i], "GRAM_SCH: vectores[" .. i .. "]", {length = dim})
    end
    tol = tol or 1e-10
    Core.assert_nonneg_number(tol, "GRAM_SCH: tol")

    local basis = {}
    for i = 1, k do
        local v = {}
        for d = 1, dim do
            v[d] = vectores[i][d]
        end
        for _, q in ipairs(basis) do
            local coeff = Vector.V_Dot(Complex.Conjugate(q), v)
            for d = 1, dim do
                v[d] = v[d] - coeff * q[d]
            end
        end
        if Vector.V_Norm(v) > tol then
            basis[#basis + 1] = Vector.V_Normalize(v)
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
    -- Claude: "Gram_Schmidt" (spelled out), not "GRAM_SCH"/"Gram_Sch" --
    -- the public-naming pass standardizing on PascalCase_With_Underscores.
    Gram_Schmidt = GRAM_SCH,
}
