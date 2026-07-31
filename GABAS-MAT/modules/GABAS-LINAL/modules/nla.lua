-- [NLA]
-- Numerical linear algebra: iterative/approximate methods, as opposed to
-- Eigen's exact/formal eigenstructure. Rayleigh moved here from eigen.lua
-- for exactly that reason -- it's an iterative numerical method, not an
-- exact decomposition.

local Core = require("GABAS-LINAL.modules.core")
local Complex = require("GABAS-LINAL.modules.complex")
local Matrix = require("GABAS-LINAL.modules.matrix")
local Vector = require("GABAS-LINAL.modules.vector")
local Systems = require("GABAS-LINAL.modules.systems")

-- Matrix-vector product Av, as a flat vector. Internal helper (used
-- by Rayleigh); not exported, since the module otherwise always represents
-- an "n x 1 matrix" as a proper matrix, not a bare vector.
local function mat_vec(A, v)
    local m, n = #A, #A[1]
    assert(#v == n, "mat_vec: vector length must match the number of columns.")
    local r = {}
    for i = 1, m do
        local row = A[i]
        local sum = 0
        for j = 1, n do
            sum = sum + row[j] * v[j]
        end
        r[i] = sum
    end
    return r
end

-- Rayleigh quotient iteration: refines a user-supplied seed vector into an
-- eigenvalue/eigenvector estimate. Starting from mu_0 = Rayleigh(v0), each
-- step solves (A - mu*I) v_new = v for a new direction and re-normalizes;
-- this converges locally cubically fast once close enough to a genuine
-- eigenpair, but it converges to WHICHEVER eigenvalue the seed vector
-- happens to be closest to -- a different seed can land on a different
-- eigenvalue, and a badly-chosen seed can (rarely) fail to converge within
-- max_iter. Works on Complex-entry matrices/seed vectors too, via the
-- Hermitian inner product (Conjugate(v) . w) used for the quotient itself --
-- unlike Eigenvalues/Eigenvectors, this never touches qr_decompose.
local function Rayleigh(matriz, v0, max_iter, tol)
    local n = #matriz
    assert(n == #matriz[1], "Rayleigh: matrix must be square.")
    assert(type(v0) == "table" and #v0 == n,
        "Rayleigh: seed vector length must match the matrix dimension.")
    assert(Vector.VNorm(v0) > Core.EPSILON, "Rayleigh: seed vector must not be the zero vector.")
    max_iter = max_iter or 100
    assert(type(max_iter) == "number" and max_iter >= 1 and max_iter == math.floor(max_iter),
        "Rayleigh: max_iter must be a positive integer.")
    tol = tol or 1e-12

    local v = Vector.VNormalize(v0)
    local mu = Vector.Vdot(Complex.Conjugate(v), mat_vec(matriz, v))
    local converged = false

    for _ = 1, max_iter do
        local shifted = Matrix.Mat_sub(matriz, Matrix.Scalar_mul(Matrix.Eye(n), mu))
        local ok, v_next = pcall(Systems.Solve, shifted, v)
        if not ok then
            -- (A - mu*I) is numerically singular: mu is already essentially
            -- an eigenvalue -- exactly what convergence looks like, so treat
            -- it as success rather than propagating Solve's error.
            converged = true
            break
        end
        v = Vector.VNormalize(v_next)
        local mu_next = Vector.Vdot(Complex.Conjugate(v), mat_vec(matriz, v))
        local delta = mu_next - mu
        mu = mu_next
        if Core.cabs(delta) < tol then
            converged = true
            break
        end
    end

    if not converged then
        io.stderr:write(string.format(
            "Rayleigh: did not converge within %d iterations (last estimate: %s); " ..
            "try a different seed vector or a larger max_iter.\n", max_iter, tostring(mu)))
    end

    return mu, v, converged
end

-- ===========================================================================
-- SVD: full singular value decomposition, A = U*Sigma*V^H, via the classical
-- Golub-Kahan-Reinsch algorithm -- Householder bidiagonalization followed by
-- implicit-shift QR iteration on the (real) bidiagonal core. U is m x m,
-- V is n x n, Sigma is m x n (both FULL/complete, not the economy/thin
-- size), all singular values nonnegative and sorted descending.
--
-- SCOPE DECISION, stated up front rather than silently: the spec this was
-- built from prefers divide-and-conquer for the bidiagonal-SVD phase, with
-- implicit QR iteration as an explicitly-sanctioned fallback ("a serious
-- implementation should value reliability over ideology"). Divide-and-
-- conquer's merge step requires safeguarded root-finding on a secular
-- equation -- a materially larger, more failure-prone undertaking than
-- everything else in this module, and not something to hand-roll without
-- a reference implementation to validate against. This implementation is
-- the QR-iteration path: a real, complete, independently-verified Full SVD
-- (every claim below was checked against known transforms and reconstructed
-- residuals with a real Lua interpreter, not just asserted), just not the
-- asymptotically-faster-for-huge-n variant.
--
-- KNOWN LIMITATION, likewise stated rather than hidden: an EXACT zero
-- appearing mid-bidiagonal (not just a small singular value, an exact zero
-- diagonal entry with nonzero neighbors) is a documented hard case for
-- plain implicit-shift QR -- real implementations (LAPACK's dbdsqr) use a
-- dedicated "zero-shift chase" for it, which is out of scope here. This
-- implementation instead detects non-convergence via max_iter and reports
-- it as an error rather than silently returning a wrong answer.
-- ===========================================================================

local function is_complex(x) return Complex.Is_complex(x) end
local function scalar_conj(x) return is_complex(x) and x:conj() or x end
local function scalar_re(x) return is_complex(x) and x.re or x end

local function conj_transpose(M)
    return Complex.Conjugate(Matrix.T(M))
end

-- Given a vector x (length k >= 1, real or Complex entries), builds a
-- Householder reflector H = I - tau*v*v^H (v[1]=1) such that H*x = beta*e1
-- with beta REAL (nonneg or negative -- sign is corrected later in
-- postprocessing) -- what lets bidiagonalizing a COMPLEX matrix produce a
-- REAL bidiagonal core, exactly what Sigma requires. Derivation: with
-- alpha=x[1], r=||x||, beta=-sign(Re(alpha))*r, v_raw=x-beta*e1,
-- w=v_raw^H*x, tau=1/w (unnormalized), then v_raw is rescaled so v[1]=1
-- and tau rescaled by |v_raw[1]|^2 to compensate -- verified directly
-- (H*x=beta*e1 and H^H*H=I, both to machine precision) against real and
-- Complex test vectors before ever being wired into bidiagonalization.
local function householder(x)
    local k = #x
    local tailnormsq = 0
    for i = 2, k do tailnormsq = tailnormsq + Core.cabs(x[i])^2 end
    local alpha = x[1]
    local r = math.sqrt(Core.cabs(alpha)^2 + tailnormsq)
    if r < Core.EPSILON then
        local v = {}
        for i = 1, k do v[i] = 0 end
        v[1] = 1
        return v, 0, 0
    end
    local sign = (scalar_re(alpha) >= 0) and 1 or -1
    local beta = -sign * r
    local v_raw = {}
    for i = 1, k do v_raw[i] = x[i] end
    v_raw[1] = v_raw[1] - beta
    local c = v_raw[1]
    local w = 0
    for i = 1, k do w = w + scalar_conj(v_raw[i]) * x[i] end
    local tau0 = 1 / w
    local v = {}
    for i = 1, k do v[i] = v_raw[i] / c end
    local tau = tau0 * (Core.cabs(c)^2)
    return v, tau, beta
end

-- Embeds (I - tau*v*v^H) (size #v) into an otherwise-identity `size x size`
-- matrix, occupying rows/cols [offset+1 .. offset+#v].
local function embed_householder(v, tau, size, offset)
    local H = Matrix.Eye(size)
    local k = #v
    for i = 1, k do
        for j = 1, k do
            local id = (i == j) and 1 or 0
            H[offset + i][offset + j] = id - tau * v[i] * scalar_conj(v[j])
        end
    end
    return H
end

-- Reduces A (m x n, m >= n) to real upper-bidiagonal B (same shape) via
-- alternating Householder reflectors, returning B, U1 (m x m), V1 (n x n)
-- with A = U1*B*V1^H. Requires m >= n (SVD() below transposes first when
-- it isn't, so this never has to handle the wide case itself).
local function bidiagonalize(A)
    local m, n = #A, #A[1]
    local Acur = Core.copy_matrix(A)
    local U1 = Matrix.Eye(m)
    local V1 = Matrix.Eye(n)

    for k = 1, n do
        -- Even a length-1 tail (x has exactly one entry, at the very last
        -- column when m==n) still needs its own reflector: for real input
        -- it's just a sign flip, but for Complex input it's the only thing
        -- forcing that last entry's PHASE to be real, and nothing later
        -- touches it to fix that up otherwise. Skipping length-1 tails was
        -- an earlier bug here -- confirmed by real random-complex-matrix
        -- testing, not just reasoned about: without this, the last diagonal
        -- (or last superdiagonal, see below) entry came out genuinely
        -- complex, not just real-with-floating-point-noise.
        local x = {}
        for i = k, m do x[#x + 1] = Acur[i][k] end
        local v, tau = householder(x)
        if tau ~= 0 then
            local H = embed_householder(v, tau, m, k - 1)
            Acur = Matrix.Mat_mul(H, Acur)
            U1 = Matrix.Mat_mul(U1, conj_transpose(H))
        end
        if k <= n - 1 then
            local y = {}
            for j = k + 1, n do y[#y + 1] = Acur[k][j] end
            local v2, tau2 = householder(y)
            if tau2 ~= 0 then
                -- Right-side reflector satisfies y_row*H_R = beta*e1_row,
                -- i.e. H_R^T (plain transpose, not conjugate) * y_col =
                -- beta*e1; algebra gives H_R = H^T = I - tau*conj(v)*conj(v)^H
                -- -- same tau, v elementwise conjugated (verified below,
                -- same way as the reflector itself was verified).
                local vbar = {}
                for i = 1, #v2 do vbar[i] = scalar_conj(v2[i]) end
                local Hr = embed_householder(vbar, tau2, n, k)
                Acur = Matrix.Mat_mul(Acur, Hr)
                V1 = Matrix.Mat_mul(V1, Hr)
            end
        end
    end
    return Acur, U1, V1
end

-- [c s; -s c] * [a;b] = [r;0], r = sqrt(a^2+b^2). LEFT-multiplication
-- convention on a column vector; a`embed_rot` built from this must be
-- transposed (negate s) before it's valid for a RIGHT multiplication --
-- this distinction is exactly what the bidiagonal QR sweep below got wrong
-- on the first attempt (verified by direct reconstruction against a known
-- eigenvalue-based cross-check before trusting it).
local function real_givens(a, b)
    if b == 0 then return 1, 0 end
    if a == 0 then return 0, 1 end
    local r
    if math.abs(a) > math.abs(b) then
        local t = b / a
        r = math.abs(a) * math.sqrt(1 + t * t)
    else
        local t = a / b
        r = math.abs(b) * math.sqrt(1 + t * t)
    end
    return a / r, b / r
end

local function embed_real_rot(c, s, size, i, j)
    local G = Matrix.Eye(size)
    G[i][i] = c;  G[i][j] = s
    G[j][i] = -s; G[j][j] = c
    return G
end

-- One Golub-Kahan implicit-shift QR step over the active block [lo,hi] of
-- a REAL n x n bidiagonal (diag d, superdiag e), n = full matrix size
-- (Ub/Vb always stay n x n; [lo,hi] just bounds which columns/rows this
-- particular sweep is allowed to touch). Wilkinson shift from the trailing
-- 2x2 of T=B^T*B; the (y,z) pair drives a right rotation (zeroing into a
-- subdiagonal bulge), then a left rotation eliminates that bulge (creating
-- the next bulge two columns over, feeding the next k). d, e mutated in
-- place; returns updated Ub, Vb.
local function golub_kahan_step(d, e, n, lo, hi, Ub, Vb)
    local dn, dn1 = d[hi], d[hi - 1]
    local en1 = e[hi - 1]
    local en2 = (hi - 2 >= lo) and e[hi - 2] or 0
    local Tnn = dn * dn + en1 * en1
    local Tn1n1 = dn1 * dn1 + en2 * en2
    local Tn1n = dn1 * en1
    local dd = (Tn1n1 - Tnn) / 2
    local mu
    if dd == 0 and Tn1n == 0 then
        mu = Tnn
    else
        local sign = (dd >= 0) and 1 or -1
        mu = Tnn - (Tn1n * Tn1n) / (dd + sign * math.sqrt(dd * dd + Tn1n * Tn1n))
    end

    local Bwork = Matrix.Zeroes(n, n)
    for i = 1, n do Bwork[i][i] = d[i] end
    for i = 1, n - 1 do Bwork[i][i + 1] = e[i] end

    local y = d[lo] * d[lo] - mu
    local z = d[lo] * e[lo]

    for k = lo, hi - 1 do
        local c, s = real_givens(y, z)
        local G1 = embed_real_rot(c, -s, n, k, k + 1)   -- transposed: right-mult
        Bwork = Matrix.Mat_mul(Bwork, G1)
        Vb = Matrix.Mat_mul(Vb, G1)

        local c2, s2 = real_givens(Bwork[k][k], Bwork[k + 1][k])
        local G2 = embed_real_rot(c2, s2, n, k, k + 1)  -- as-is: left-mult
        Bwork = Matrix.Mat_mul(G2, Bwork)
        Ub = Matrix.Mat_mul(Ub, Matrix.T(G2))

        if k < hi - 1 then
            y = Bwork[k][k + 1]
            z = Bwork[k][k + 2]
        end
    end

    for i = 1, n do d[i] = Bwork[i][i] end
    for i = 1, n - 1 do e[i] = Bwork[i][i + 1] end
    return Ub, Vb
end

-- Diagonalizes a real n x n bidiagonal (d,e) via repeated Golub-Kahan
-- steps with Demmel-Kahan relative deflation (|e_i| <= tol*(|d_i|+|d_{i+1}|)
-- splits the problem). Returns the (still signed, unsorted) diagonal and
-- the accumulated real orthogonal Ub, Vb with diag(d) = Ub^T*B*Vb... i.e.
-- B = Ub*diag(d)*Vb^T. Raises an error (does not silently return a wrong
-- answer) if max_iter sweeps pass without full deflation -- see the
-- "KNOWN LIMITATION" note on SVD above for when this can happen.
local function bidiagonal_svd(d_in, e_in, max_iter)
    local n = #d_in
    local d, e = {}, {}
    for i = 1, n do d[i] = d_in[i] end
    for i = 1, n - 1 do e[i] = e_in[i] end
    local Ub, Vb = Matrix.Eye(n), Matrix.Eye(n)
    if n == 1 then return d, Ub, Vb end

    local hi = n
    local sweeps = 0
    while hi > 1 do
        sweeps = sweeps + 1
        assert(sweeps <= max_iter,
            "SVD: bidiagonal QR iteration did not converge within max_iter sweeps " ..
            "(possible exact-zero-mid-diagonal case -- see the module's documented " ..
            "limitation; try increasing max_iter, or note that the matrix may need " ..
            "a dedicated rank-revealing approach this implementation doesn't cover).")
        local lo = hi
        while lo > 1 do
            local tol = Core.EPSILON * (math.abs(d[lo - 1]) + math.abs(d[lo]))
            if math.abs(e[lo - 1]) <= tol then
                e[lo - 1] = 0
                break
            end
            lo = lo - 1
        end
        if lo == hi then
            hi = hi - 1
        else
            Ub, Vb = golub_kahan_step(d, e, n, lo, hi, Ub, Vb)
            local tol2 = Core.EPSILON * (math.abs(d[hi - 1]) + math.abs(d[hi]))
            if math.abs(e[hi - 1]) <= tol2 then
                e[hi - 1] = 0
                hi = hi - 1
            end
        end
    end
    return d, Ub, Vb
end

-- Full Singular Value Decomposition: A = U * Sigma * V^H.
--   A: m x n matrix (real numbers, Complex values, or a mix).
--   opts (optional table):
--     max_iter (default 500): bidiagonal QR sweep cap, see bidiagonal_svd.
--     verify (default false): if true, additionally asserts the residual
--       ||A-U*Sigma*V^H||, ||U^H U - I||, ||V^H V - I|| are all small and
--       that singular values are sorted descending nonneg -- the spec's
--       own optional Phase 7/Step 10 verification, exposed as an opt-in
--       rather than paid on every call.
-- Returns U (m x m), Sigma (m x n), V (n x n), sv (flat array of the
-- min(m,n) singular values, already sorted descending -- the same
-- Sigma[i][i] values, handed back the way Eigenvalues() hands back
-- eigenvalues, so comparing/thresholding them doesn't require digging
-- them out of Sigma's diagonal yourself).
local function SVD(A, opts)
    assert(type(A) == "table" and #A > 0 and type(A[1]) == "table" and #A[1] > 0,
        "SVD: A must be a non-empty matrix (a table of non-empty row-tables).")
    local m, n = #A, #A[1]
    for i = 1, m do
        assert(type(A[i]) == "table" and #A[i] == n,
            "SVD: A must be rectangular -- every row must have the same length.")
        for j = 1, n do
            local v = A[i][j]
            if is_complex(v) then
                assert(type(v.re) == "number" and type(v.im) == "number" and
                    v.re == v.re and v.im == v.im and
                    v.re ~= math.huge and v.re ~= -math.huge and
                    v.im ~= math.huge and v.im ~= -math.huge,
                    "SVD: entry A[" .. i .. "][" .. j .. "] is not finite (NaN/Inf).")
            else
                assert(type(v) == "number",
                    "SVD: entry A[" .. i .. "][" .. j .. "] must be a number or a Complex value.")
                assert(v == v and v ~= math.huge and v ~= -math.huge,
                    "SVD: entry A[" .. i .. "][" .. j .. "] is not finite (NaN/Inf).")
            end
        end
    end
    opts = opts or {}
    local max_iter = opts.max_iter or 500
    assert(type(max_iter) == "number" and max_iter >= 1 and max_iter == math.floor(max_iter),
        "SVD: opts.max_iter must be a positive integer.")

    -- All-zero matrix: quick exit, per the spec's Phase 2 -- also sidesteps
    -- a division by zero in the scaling step below.
    local alpha = 0
    for i = 1, m do
        for j = 1, n do
            local v = Core.cabs(A[i][j])
            if v > alpha then alpha = v end
        end
    end
    if alpha == 0 then
        return Matrix.Eye(m), Matrix.Zeroes(m, n), Matrix.Eye(n)
    end

    -- Safe scaling: brings the largest-magnitude entry to 1 before doing
    -- any arithmetic, so Householder norms/rotations stay well away from
    -- overflow/underflow regardless of A's own scale; unscaled at the end.
    local SAFE_MIN, SAFE_MAX = 1e-100, 1e100
    local scale = 1
    if alpha < SAFE_MIN or alpha > SAFE_MAX then
        scale = 1 / alpha
    end
    local Ascaled = (scale == 1) and A or Matrix.Scalar_mul(A, scale)

    -- Core algorithm assumes "tall" (rows >= cols); for a wide A, solve
    -- for A^H instead (now tall) and un-transpose the result: if
    -- A^H = U'*Sigma'*V'^H then A = (A^H)^H = V'*Sigma'^T*U'^H, i.e.
    -- U=V', Sigma=Sigma'^T (just the transposed shape), V=U'.
    local transposed = m < n
    local Awork = transposed and conj_transpose(Ascaled) or Ascaled
    local mw, nw = #Awork, #Awork[1]

    local B, U1, V1 = bidiagonalize(Awork)
    -- B is mathematically real by construction (the Householder reflectors
    -- were built specifically to make it so, see householder() above), but
    -- Complex arithmetic never auto-demotes to a plain number even when the
    -- imaginary part is exactly 0 -- so a chain of Complex arithmetic can
    -- still hand back e.g. Complex(5, 0) instead of plain 5. Extract the
    -- real part explicitly here (bidiagonal_svd below is real-arithmetic
    -- only); assert the imaginary part actually IS negligible first, since
    -- if it weren't, that would mean an actual bug upstream, not just a
    -- representation quirk.
    local function assert_real(x, where)
        if is_complex(x) then
            assert(Core.cabs(x.im) < 1e-6 * (Core.cabs(x.re) + 1),
                "SVD: internal error -- expected a real bidiagonal entry (" ..
                where .. ") but got " .. tostring(x) .. ".")
            return x.re
        end
        return x
    end
    local d, e = {}, {}
    for i = 1, nw do d[i] = assert_real(B[i][i], "diagonal") end
    for i = 1, nw - 1 do e[i] = assert_real(B[i][i + 1], "superdiagonal") end

    local dsvd, Ub, Vb = bidiagonal_svd(d, e, max_iter)

    -- Force singular values nonnegative: flip the sign onto the
    -- corresponding column of Ub instead.
    for i = 1, nw do
        if dsvd[i] < 0 then
            dsvd[i] = -dsvd[i]
            for r = 1, nw do Ub[r][i] = -Ub[r][i] end
        end
    end

    -- Sort descending, permuting Ub/Vb columns to match.
    local order = {}
    for i = 1, nw do order[i] = i end
    table.sort(order, function(a, b) return dsvd[a] > dsvd[b] end)
    local dsorted, Ubsorted, Vbsorted = {}, Matrix.Zeroes(nw, nw), Matrix.Zeroes(nw, nw)
    for newi, oldi in ipairs(order) do
        dsorted[newi] = dsvd[oldi]
        for r = 1, nw do
            Ubsorted[r][newi] = Ub[r][oldi]
            Vbsorted[r][newi] = Vb[r][oldi]
        end
    end

    -- Embed Ub (nw x nw) into the top-left of an mw x mw identity: rows
    -- nw+1..mw of B are already all-zero, so any transformation restricted
    -- to Ub's own block leaves them untouched.
    local UbFull = Matrix.Eye(mw)
    for i = 1, nw do
        for j = 1, nw do UbFull[i][j] = Ubsorted[i][j] end
    end

    local U = Matrix.Mat_mul(U1, UbFull)
    local V = Matrix.Mat_mul(V1, Vbsorted)
    local Sigma = Matrix.Zeroes(mw, nw)
    for i = 1, nw do Sigma[i][i] = dsorted[i] / scale end

    -- Flat, already-sorted-descending singular values -- min(m,n) of them,
    -- the same convention Eigenvalues() already uses (a plain array, not
    -- values buried on a matrix diagonal you'd have to dig out yourself).
    -- Unaffected by the tall/wide transpose trick below: the singular
    -- values themselves don't change, only which of U/V they came from.
    local sv = {}
    for i = 1, nw do sv[i] = dsorted[i] / scale end

    if transposed then
        U, V = V, U
        Sigma = Matrix.T(Sigma)
    end

    if opts.verify then
        local recon = Matrix.Mat_mul(Matrix.Mat_mul(U, Sigma), conj_transpose(V))
        local resid = 0
        for i = 1, m do
            for j = 1, n do
                local diff = Core.cabs(recon[i][j] - A[i][j])
                if diff > resid then resid = diff end
            end
        end
        local tol = 1e-8 * (alpha + 1)
        assert(resid < tol,
            "SVD: verification failed -- ||A - U*Sigma*V^H||_max = " .. resid ..
            " exceeds tolerance " .. tol .. ".")
        local function orth_err(M)
            local n2 = #M[1]
            local MhM = Matrix.Mat_mul(conj_transpose(M), M)
            local worst = 0
            for i = 1, n2 do
                for j = 1, n2 do
                    local id = (i == j) and 1 or 0
                    local diff = Core.cabs(MhM[i][j] - id)
                    if diff > worst then worst = diff end
                end
            end
            return worst
        end
        assert(orth_err(U) < 1e-8, "SVD: verification failed -- U is not unitary.")
        assert(orth_err(V) < 1e-8, "SVD: verification failed -- V is not unitary.")
        local p = math.min(m, n)
        for i = 1, p - 1 do
            assert(Sigma[i][i] >= Sigma[i + 1][i + 1] - 1e-12,
                "SVD: verification failed -- singular values are not sorted descending.")
        end
        assert(p == 0 or Sigma[p][p] >= -1e-12,
            "SVD: verification failed -- a singular value came out negative.")
    end

    return U, Sigma, V, sv
end

return {
    Rayleigh = Rayleigh,
    SVD = SVD,
}

-- B"H.