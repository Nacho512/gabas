function Show_matrix(matriz)
    local n = #matriz[1]
    local m = #matriz
    for i = 1, m do
        for j = 1, n do
            io.write(matriz[i][j], "\t")
        end
        print("\n")
    end
end

function Matrix(data, num_filas)
    local limite = #data
    assert((limite % num_filas)==0, "Mismatch between number of elements and number of rows.")
    local m = {}
    local prop = limite/num_filas
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

function Sequence(aleph,tat)
    local m ={}
    local k = 1
    for j = aleph, tat do
        m[k] = j
        k = k + 1
    end
    return m
end

function Show_table(tabla)
    local limite = #tabla
    for j = 1, limite do
        print(tabla[j])
    end
end

function Mat_mul(M1,M2)
    local m = #M1
    local p = #M1[1]
    local n = #M2[1]
    local C = {}
    for i=1, m do
        C[i] = {}
        for j=1, n do
            local sum = 0
            for k =1, p do
                sum = sum + M1[i][k]*M2[k][j]
            end
            C[i][j] = sum
        end
    end
    return C
end

function Mat_sum(M1,M2)
    local m = #M1
    local n = #M1[1]
    local C = {}
    for i=1, m do
        C[i] = {}
        for j = 1, n do
            C[i][j] = M1[i][j] + M2[i][j]
        end
    end
    return C
end

function Eye(dim)
    local C = {}
    for i = 1, dim do
        C[i] = {}
        for j =1, dim do
            if i == j then
                C[i][j] = 1
            else
                C[i][j] = 0
            end
        end
    end
    return C
end

function Zeroes(n_rows, n_cols)
    local C = {}
    for i = 1, n_rows do
        C[i] = {}
        for j =1, n_cols do
            C[i][j] = 0
        end
    end
    return C
end

function Random_mat(m,n,valor)
    local C = {}
    for i = 1, m do
        C[i] = {}
        for j =1, n do
            C[i][j] = math.random(1,valor)
        end
    end
    return C
end

function T(matriz)
    local m = #matriz
    local n = #matriz[1]
    local C = {}
    for i = 1, n do
        C[i] = {}
        for j = 1, m do
            C[i][j] = matriz[j][i]
        end
    end
    return C
end

function Determinant(matrix)
    local rows, cols = #matrix, #matrix[1]
    assert(rows == cols, "Mismatch between number of elements and number of rows.")
    
    -- LU decomposition
    local lu = {}
    for i = 1, rows do
        lu[i] = {}
        for j = 1, cols do
            lu[i][j] = matrix[i][j]
        end
    end
    
    for k = 1, rows - 1 do
        for i = k + 1, rows do
            local factor = lu[i][k] / lu[k][k]
            for j = k + 1, rows do
                lu[i][j] = lu[i][j] - factor * lu[k][j]
            end
            lu[i][k] = factor
        end
    end
    
    -- determinant calculation
    local det = 1
    for i = 1, rows do
        det = det * lu[i][i]
    end
    
    return det
end

return {Show_matrix=Show_matrix, Matrix=Matrix, Sequence=Sequence, Show_table=Show_table, Mat_mul=Mat_mul, Mat_sum=Mat_sum,
Eye=Eye, Zeroes=Zeroes, Random_mat=Random_mat, T=T, Determinant=Determinant}