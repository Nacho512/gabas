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

local function Determinant(matrix)
    local rows, cols = #matrix, #matrix[1]
    assert(rows == cols, "Determinant: matrix must be square.")

    -- LU decomposition (Gaussian elimination) with partial pivoting.
    local lu = {}
    for i = 1, rows do
        lu[i] = {}
        for j = 1, cols do
            lu[i][j] = matrix[i][j]
        end
    end

    local sign = 1
    for k = 1, rows - 1 do
        if lu[k][k] == 0 then
            local pivot_row = nil
            for i = k + 1, rows do
                if lu[i][k] ~= 0 then
                    pivot_row = i
                    break
                end
            end
            if not pivot_row then
                -- Entire column below (and at) the pivot is zero: matrix is singular.
                return 0
            end
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

    local det = sign
    for i = 1, rows do
        det = det * lu[i][i]
    end

    return det
end

return {Show_matrix=Show_matrix, Matrix=Matrix, Sequence=Sequence, Show_table=Show_table, Mat_mul=Mat_mul, Mat_sum=Mat_sum,
Eye=Eye, Zeroes=Zeroes, Random_mat=Random_mat, T=T, Determinant=Determinant}
