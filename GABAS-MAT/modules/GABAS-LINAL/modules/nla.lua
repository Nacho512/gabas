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

return {
    Rayleigh = Rayleigh,
}

-- B"H.