-- [Matrix]
-- Construction, display, and arithmetic for 2-D matrices. Self-contained
-- by design: none of these functions call Is_complex/cabs/etc. directly --
-- Complex-entry support flows through automatically via the +/-/* operator
-- overloads defined in Complex (see that module's header comment), so this
-- module doesn't need to require Complex at all.

local function Show_matrix(matriz)
    assert(type(matriz) == "table" and #matriz > 0 and type(matriz[1]) == "table" and #matriz[1] > 0,
        "Show_matrix: matriz must be a non-empty matrix (a table of non-empty row-tables).")
    local m, n = #matriz, #matriz[1]
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
    assert(type(num_filas) == "number" and num_filas >= 1 and num_filas == math.floor(num_filas),
        "Matrix: num_filas must be a positive integer.")
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
    assert(type(aleph) == "number" and type(tat) == "number", "Sequence: aleph and tat must be numbers.")
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
    assert(type(dim) == "number" and dim >= 1 and dim == math.floor(dim), "Eye: dim must be a positive integer.")
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
    assert(type(n_rows) == "number" and n_rows >= 1 and n_rows == math.floor(n_rows),
        "Zeroes: n_rows must be a positive integer.")
    assert(type(n_cols) == "number" and n_cols >= 1 and n_cols == math.floor(n_cols),
        "Zeroes: n_cols must be a positive integer.")
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
    assert(type(m) == "number" and m >= 1 and m == math.floor(m), "Random_mat: m must be a positive integer.")
    assert(type(n) == "number" and n >= 1 and n == math.floor(n), "Random_mat: n must be a positive integer.")
    assert(type(valor) == "number" and valor >= 1, "Random_mat: valor must be >= 1.")
    if seed ~= nil then
        assert(type(seed) == "number", "Random_mat: seed must be a number.")
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
    assert(type(matriz) == "table" and #matriz > 0 and type(matriz[1]) == "table" and #matriz[1] > 0,
        "T: matriz must be a non-empty matrix.")
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

return {
    Show_matrix = Show_matrix,
    Matrix = Matrix,
    Sequence = Sequence,
    Show_table = Show_table,
    Mat_mul = Mat_mul,
    Mat_sum = Mat_sum,
    Mat_sub = Mat_sub,
    Scalar_mul = Scalar_mul,
    Eye = Eye,
    Zeroes = Zeroes,
    Random_mat = Random_mat,
    T = T,
}
