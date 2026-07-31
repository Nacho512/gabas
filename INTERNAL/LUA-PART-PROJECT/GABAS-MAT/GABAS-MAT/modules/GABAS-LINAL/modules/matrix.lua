-- [Matrix]
-- Construction, display, and arithmetic for 2-D matrices.
-- Codex: Core supplies the shared public-input validation policy; complex
-- arithmetic still flows through the operators defined by Complex.
local Core = require("GABAS-LINAL.modules.core")
local Complex = require("GABAS-LINAL.modules.complex")

local function Show_matrix(matriz)
    local m, n = Core.assert_matrix(matriz, "Show_matrix: matriz")
    for i = 1, m do
        local row = matriz[i]
        for j = 1, n do
            -- tostring() (not io.write's own number coercion) so Complex
            -- entries render via __tostring; io.write itself only accepts
            -- plain strings/numbers and would error on a Complex table.
            io.write(tostring(row[j]), "\t")
        end
        print()
    end
end

local function Matrix(data, num_filas)
    assert(type(data) == "table" and #data > 0, "Matrix: data must be a non-empty table.")
    Core.assert_positive_integer(num_filas, "Matrix: num_filas")
    local limite = #data
    assert((limite % num_filas) == 0, "Matrix: mismatch between number of elements and number of rows.")
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
    Core.assert_finite_number(aleph, "Sequence: aleph")
    Core.assert_finite_number(tat, "Sequence: tat")
    local m = {}
    local k = 1
    for j = aleph, tat do
        m[k] = j
        k = k + 1
    end
    return m
end

local function Show_table(tabla)
    assert(type(tabla) == "table", "Show_table: tabla must be a table.")
    local limite = #tabla
    for j = 1, limite do
        print(tabla[j])
    end
end

local function Mat_mul(M1, M2)
    local m, p = Core.assert_matrix(M1, "Mat_mul: M1")
    local m2, n = Core.assert_matrix(M2, "Mat_mul: M2")
    assert(m2 == p, "Mat_mul: number of columns in M1 must match number of rows in M2.")
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
    local m, n = Core.assert_matrix(M1, "Mat_sum: M1")
    local m2, n2 = Core.assert_matrix(M2, "Mat_sum: M2")
    assert(m2 == m and n2 == n, "Mat_sum: matrices must have matching dimensions.")
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
    local m, n = Core.assert_matrix(M1, "Mat_sub: M1")
    local m2, n2 = Core.assert_matrix(M2, "Mat_sub: M2")
    assert(m2 == m and n2 == n, "Mat_sub: matrices must have matching dimensions.")
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
    local m, n = Core.assert_matrix(matriz, "Scalar_mul: matriz")
    Core.assert_finite_scalar(escalar, "Scalar_mul: escalar")
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
    Core.assert_positive_integer(dim, "Eye: dim")
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
    Core.assert_positive_integer(n_rows, "Zeroes: n_rows")
    Core.assert_positive_integer(n_cols, "Zeroes: n_cols")
    local C = {}
    for i = 1, n_rows do
        C[i] = {}
        for j = 1, n_cols do
            C[i][j] = 0
        end
    end
    return C
end

-- If `seed` is given, reseeds Lua's global PRNG right before filling the
-- matrix, so the same seed always reproduces the same matrix -- useful for
-- reproducible tests and examples, and a version-independent fix for the
-- fact that math.random's underlying generator (and its default seeding
-- quality) differs across Lua versions: pre-5.4 delegates straight to the
-- C library's rand() with "no guarantees on its statistical properties"
-- (the Lua 5.1 manual's own words), while 5.4 switched to xoshiro256**.
-- Omitting `seed` leaves the current PRNG state exactly as it is, rather
-- than reseeding on every call -- reseeding unconditionally here would be
-- worse than doing nothing, since two calls made in the same run (or two
-- runs started in the same second, on a Lua where the default seed comes
-- from time()) could otherwise reset to the same state and produce
-- identical "random" matrices.
local function Random_mat(m, n, valor, seed)
    Core.assert_positive_integer(m, "Random_mat: m")
    Core.assert_positive_integer(n, "Random_mat: n")
    Core.assert_positive_integer(valor, "Random_mat: valor")
    if seed ~= nil then
        Core.assert_integer(seed, "Random_mat: seed")
        math.randomseed(seed)
    end
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
    local m, n = Core.assert_matrix(matriz, "T: matriz")
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

local function Matrix_Vector_Mul(matriz, vector)
    local rows, cols = Core.assert_matrix(matriz, "Matrix_Vector_Mul: matriz")
    Core.assert_vector(vector, "Matrix_Vector_Mul: vector", {length = cols})
    local result = {}
    for i = 1, rows do
        local sum = 0
        for j = 1, cols do sum = sum + matriz[i][j] * vector[j] end
        result[i] = sum
    end
    return result
end

local function Conjugate_Transpose(matriz)
    Core.assert_matrix(matriz, "Conjugate_Transpose: matriz")
    return Complex.Conjugate(T(matriz))
end

return {
    Show_matrix = Show_matrix,
    Show_Matrix = Show_matrix,
    Matrix = Matrix,
    Sequence = Sequence,
    Show_table = Show_table,
    Show_Table = Show_table,
    Mat_mul = Mat_mul,
    Mat_Mul = Mat_mul,
    Mat_sum = Mat_sum,
    Mat_Sum = Mat_sum,
    Mat_sub = Mat_sub,
    Mat_Sub = Mat_sub,
    Scalar_mul = Scalar_mul,
    Scalar_Mul = Scalar_mul,
    Eye = Eye,
    Zeroes = Zeroes,
    Random_mat = Random_mat,
    Random_Mat = Random_mat,
    T = T,
    Matrix_Vector_Mul = Matrix_Vector_Mul,
    Mat_vec = Matrix_Vector_Mul,
    Conjugate_Transpose = Conjugate_Transpose,
}
