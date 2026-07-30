-- ===========================================================================
-- SUBMODULOS TENTATIVOS
-- Cada funcion de este archivo lleva una etiqueta [Categoria] al inicio de
-- su comentario, anticipando a que submodulo se mudaria el dia que este
-- archivo se separe en varios. Es solo documentacion por ahora -- no cambia
-- como se cargan ni como se llaman las funciones todavia.
--
--   Core        -- infraestructura interna compartida por 2+ categorias
--                  (EPSILON, cabs, copy_matrix, find_pivot_row, next_pow2)
--   Complex     -- Complex, sus metametodos/metodos, Conjugate
--   Matrix      -- construccion/despliegue/aritmetica de matrices 2-D
--   Tensor      -- generalizacion n-dimensional de la familia Matrix
--   Vector      -- VTrace, Vdot, VNorm, VNormalize
--   Determinant -- MDet, Big_MDet y el multiplicador Strassen que usa
--   Systems     -- Rank, Inverse, Solve
--   Eigen       -- Eigenvalues, Eigenvectors, Rayleigh, GRAM_SCH
--   FFT         -- FFT, IFFT
-- ===========================================================================

-- [Core]
local EPSILON = 1e-9

-- ===========================================================================
-- COMPLEX NUMBERS
-- A Complex value is a table {re=..., im=...} with a metatable implementing
-- +, -, *, /, unary -, ==, and tostring. Because Mat_mul/Mat_sum/Mat_sub/
-- Scalar_mul/Vdot/VTrace/MDet/Rank/Inverse/Solve are all written using
-- plain +, -, *, / on their entries (no math.* calls baked in), they gain
-- complex-number support "for free" through these operators -- Lua looks up
-- __add/__sub/__mul/__div on whichever operand is a table, and each
-- metamethod promotes a plain number operand to Complex(x, 0) automatically.
-- A matrix/vector with only plain-number entries is completely unaffected:
-- the metamethods are never invoked and ordinary Lua number arithmetic runs.
-- ===========================================================================

-- [Complex]
local Complex_mt = {}
Complex_mt.__index = Complex_mt

-- [Complex]
local function is_complex(x)
    return type(x) == "table" and getmetatable(x) == Complex_mt
end

-- [Complex]
local function Complex(re, im)
    re, im = re or 0, im or 0
    assert(type(re) == "number" and type(im) == "number", "Complex: re and im must be numbers.")
    return setmetatable({re = re, im = im}, Complex_mt)
end

-- [Complex]
local function to_complex(x)
    if is_complex(x) then return x end
    assert(type(x) == "number", "Complex: expected a number or a Complex value.")
    return Complex(x, 0)
end

-- [Complex] operator metamethods (__add.._tostring) and instance methods
-- (conj/abs/arg) for Complex_mt -- see the header comment above for why
-- these are what let real arithmetic gain Complex support "for free".
function Complex_mt.__add(a, b)
    a, b = to_complex(a), to_complex(b)
    return Complex(a.re + b.re, a.im + b.im)
end

function Complex_mt.__sub(a, b)
    a, b = to_complex(a), to_complex(b)
    return Complex(a.re - b.re, a.im - b.im)
end

function Complex_mt.__mul(a, b)
    a, b = to_complex(a), to_complex(b)
    return Complex(a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re)
end

function Complex_mt.__div(a, b)
    a, b = to_complex(a), to_complex(b)
    local denom = b.re * b.re + b.im * b.im
    assert(denom > EPSILON, "Complex: division by zero.")
    return Complex((a.re * b.re + a.im * b.im) / denom, (a.im * b.re - a.re * b.im) / denom)
end

function Complex_mt.__unm(a)
    return Complex(-a.re, -a.im)
end

function Complex_mt.__eq(a, b)
    return a.re == b.re and a.im == b.im
end

function Complex_mt.__tostring(a)
    if a.im == 0 then
        return tostring(a.re)
    end
    local sign = a.im < 0 and "-" or "+"
    return tostring(a.re) .. sign .. tostring(math.abs(a.im)) .. "i"
end

function Complex_mt.conj(a)
    return Complex(a.re, -a.im)
end

function Complex_mt.abs(a)
    return math.sqrt(a.re * a.re + a.im * a.im)
end

function Complex_mt.arg(a)
    return math.atan(a.im, a.re)
end

-- [Core] Magnitude that works on both plain numbers and Complex values --
-- used everywhere a pivoting/convergence/near-zero check needs "how big is
-- this", regardless of whether the entries are real or complex. Shared by
-- Determinant, Systems, and Eigen alike, which is what makes it Core rather
-- than belonging to any one future submodule.
local function cabs(x)
    if is_complex(x) then
        return x:abs()
    end
    return math.abs(x)
end

-- [Complex] Elementwise conjugate. Accepts either a vector (flat table) or a
-- matrix (table of row-tables); plain-number entries pass through
-- unchanged. Lives conceptually in Complex even though its only current
-- callers (GRAM_SCH, Rayleigh) are in Eigen -- that cross-submodule call is
-- normal (Eigen would `require` Complex), unlike the Core helpers above,
-- which have no natural domain-category home of their own.
local function Conjugate(x)
    local C = {}
    if type(x[1]) == "table" and not is_complex(x[1]) then
        for i = 1, #x do
            local row, row_c = x[i], {}
            for j = 1, #row do
                row_c[j] = is_complex(row[j]) and row[j]:conj() or row[j]
            end
            C[i] = row_c
        end
    else
        for i = 1, #x do
            C[i] = is_complex(x[i]) and x[i]:conj() or x[i]
        end
    end
    return C
end

-- [Eigen] Matrix-vector product Av, as a flat vector. Internal helper (used
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

-- [Core] Used by MDet/Big_MDet (Determinant), Rank (Systems), and
-- Eigenvalues/Eigenvectors (Eigen) alike, so callers never mutate the
-- matrix they passed in.
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

-- [Core] Searches rows [k, last_row] of column `col` for the entry with the
-- largest magnitude (classic partial pivoting: picking the largest
-- magnitude, not just the first non-zero one, keeps the elimination
-- numerically stable). Returns the row index, or nil if every candidate is
-- ~0 (singular column). Shared by Determinant and Systems.
local function find_pivot_row(mat, k, last_row, col)
    local best_row, best_val = nil, EPSILON
    for i = k, last_row do
        local v = cabs(mat[i][col])
        if v > best_val then
            best_row, best_val = i, v
        end
    end
    return best_row
end

-- [Tensor] Shared guard for the whole *_tensor family (Tensor, Zeroes_tensor,
-- Eye_tensor, Random_tensor): `dims` must explicitly list every dimension
-- size, with at least 3 entries -- 2 or fewer belongs to the matrix-only
-- counterpart named by `counterpart_name`, since (as established for
-- Tensor) there is no unique way to infer more than one missing dimension
-- from a single count. `fn_name` and `counterpart_name` are spliced into
-- every message so the error always names the actual caller.
local function validate_tensor_dims(dims, fn_name, counterpart_name)
    assert(type(dims) == "table" and #dims >= 3,
        fn_name .. ": dims must be a table listing every dimension size, with at least 3 entries " ..
        "(use " .. counterpart_name .. " for 2 dimensions, or a plain flat table for 1).")
    for i = 1, #dims do
        local extent = dims[i]
        assert(type(extent) == "number" and extent >= 1 and extent == math.floor(extent),
            fn_name .. ": dims[" .. i .. "] must be a positive integer.")
    end
end

-- [Tensor] Builds a tensor of shape dims[dim_idx..#dims], every entry set to
-- `value`. Shared by Zeroes_tensor (value = 0) and by Eye_tensor's
-- off-diagonal branches (also value = 0, but reached from a different
-- starting dim_idx).
local function build_filled_tensor(dims, dim_idx, value)
    local extent = dims[dim_idx]
    local t = {}
    if dim_idx == #dims then
        for i = 1, extent do
            t[i] = value
        end
        return t
    end
    for i = 1, extent do
        t[i] = build_filled_tensor(dims, dim_idx + 1, value)
    end
    return t
end

-- [Matrix]
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

-- [Tensor] Counts how many levels of table-nesting `t` has before hitting a leaf
-- (a plain number or a Complex value) -- i.e. how many dimensions a tensor
-- built by Tensor/Zeroes_tensor/Eye_tensor/Random_tensor actually has.
-- Only ever walks down through t[1], t[1][1], ...: like the rest of the
-- module, a ragged (non-rectangular) tensor is assumed not to occur.
local function tensor_depth(t)
    if type(t) ~= "table" or is_complex(t) then
        return 0
    end
    return 1 + tensor_depth(t[1])
end

-- [Tensor] Prints one (depth-2)-dimensional matrix "slice" of a tensor per call,
-- each preceded by a header naming the fixed leading indices that got us
-- there (e.g. "T[2][1] =" for a slice fixing the first two of four
-- dimensions), then recurses one dimension shallower until only 2 remain.
local function show_tensor_rec(t, depth, path)
    if depth == 2 then
        io.write("T[", table.concat(path, "]["), "] =\n")
        for i = 1, #t do
            local row = t[i]
            for j = 1, #row do
                io.write(tostring(row[j]), "\t")
            end
            print()
        end
        print()
        return
    end
    for i = 1, #t do
        path[#path + 1] = i
        show_tensor_rec(t[i], depth - 1, path)
        path[#path] = nil
    end
end

-- [Tensor] Tensor counterpart of Show_matrix: prints every 2-D slice of a tensor of
-- 3 or more dimensions, each labeled with the leading indices fixed to
-- reach it, so an n-dimensional block of numbers stays readable as a
-- sequence of ordinary matrices instead of one flat wall of digits.
local function Show_tensor(t)
    assert(type(t) == "table" and #t > 0 and type(t[1]) == "table" and #t[1] > 0,
        "Show_tensor: t must be a non-empty tensor (nested tables).")
    local depth = tensor_depth(t)
    assert(depth >= 3,
        "Show_tensor: t must be a tensor of at least 3 dimensions (use Show_matrix for 2).")
    show_tensor_rec(t, depth, {})
end

-- [Matrix]
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

-- [Tensor] Recursively fills the innermost-to-outermost nesting of a tensor, one
-- dimension per recursion level, threading `offset` through so each leaf of
-- `data` is consumed exactly once, in the same row-major order Matrix uses
-- (the LAST dimension varies fastest) -- consistent with the rest of the
-- module, where a 2-D matrix is already a table of row-tables.
local function build_tensor(data, dims, dim_idx, offset)
    local extent = dims[dim_idx]
    local t = {}
    if dim_idx == #dims then
        for i = 1, extent do
            t[i] = data[offset + i]
        end
        return t, offset + extent
    end
    for i = 1, extent do
        t[i], offset = build_tensor(data, dims, dim_idx + 1, offset)
    end
    return t, offset
end

-- [Tensor] Tensor generalizes Matrix from 2 dimensions to any number of them: `data`
-- stays a flat vector (exactly like Matrix's `data`), and `dims` is a table
-- listing EVERY dimension's size explicitly -- e.g. {2, 3, 4} for a 2x3x4
-- tensor -- rather than trying to infer all-but-one dimension from a single
-- count the way Matrix infers columns from num_filas. That inference does
-- not generalize past 2 dimensions (there is no unique way to split a
-- leftover element count across 2+ unknown dimensions), so for a tensor
-- every dimension size must be given up front. The result is nested tables
-- #dims levels deep: a 3-D tensor T is indexed T[i][j][k].
local function Tensor(data, dims)
    assert(type(data) == "table" and #data > 0, "Tensor: data must be a non-empty table.")
    validate_tensor_dims(dims, "Tensor", "Matrix")
    local total = 1
    for i = 1, #dims do
        total = total * dims[i]
    end
    assert(#data == total,
        "Tensor: mismatch between number of elements (" .. #data ..
        ") and the product of dims (" .. total .. ").")
    local t = build_tensor(data, dims, 1, 0)
    return t
end

-- [Matrix]
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

-- [Matrix]
local function Show_table(tabla)
    assert(type(tabla) == "table", "Show_table: tabla must be a table.")
    local limite = #tabla
    for j = 1, limite do
        print(tabla[j])
    end
end

-- [Matrix]
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

-- [Matrix]
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

-- [Matrix]
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

-- [Matrix]
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

-- [Matrix]
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

-- [Tensor] Builds levels dim_idx..#dims of the hyperdiagonal, given that every index
-- chosen so far exactly equals ref_index (the value the very first index
-- took) -- i.e. the path built so far still sits "on the diagonal". Once a
-- branch picks an index different from ref_index, everything below it is
-- necessarily 0, so that branch is handed off to build_filled_tensor
-- instead of continuing to track a (now pointless) reference index.
local function build_eye_tensor(dims, dim_idx, ref_index)
    local extent = dims[dim_idx]
    local t = {}
    if dim_idx == #dims then
        for i = 1, extent do
            t[i] = (i == ref_index) and 1 or 0
        end
        return t
    end
    for i = 1, extent do
        if i == ref_index then
            t[i] = build_eye_tensor(dims, dim_idx + 1, ref_index)
        else
            t[i] = build_filled_tensor(dims, dim_idx + 1, 0)
        end
    end
    return t
end

-- [Tensor] Tensor counterpart of Eye. A matrix's identity property (Eye(n) is the
-- identity element for Mat_mul) does not generalize to order-3-or-higher
-- tensors -- this module never defines a general "tensor product" for
-- Mat_mul to be an identity element of in the first place. What DOES
-- generalize directly is Eye's literal defining rule -- 1 exactly where
-- every index coincides, 0 everywhere else -- which is the standard
-- "hyperdiagonal" (or super-diagonal) tensor from tensor-algebra
-- literature. Because "every index coincides" only has one obvious meaning
-- when every axis has the same length, dims must describe a hypercube
-- (all entries equal); for an irregular shape, build your own fill logic
-- on top of Zeroes_tensor instead.
local function Eye_tensor(dims)
    validate_tensor_dims(dims, "Eye_tensor", "Eye")
    for i = 2, #dims do
        assert(dims[i] == dims[1],
            "Eye_tensor: all dimensions must be equal (got dims[" .. i .. "] = " .. dims[i] ..
            ", dims[1] = " .. dims[1] .. ") -- the hyperdiagonal is only well-defined for a " ..
            "hypercubic shape; use Zeroes_tensor plus your own fill logic for an irregular one.")
    end
    local extent = dims[1]
    local t = {}
    for i = 1, extent do
        t[i] = build_eye_tensor(dims, 2, i)
    end
    return t
end

-- [Matrix]
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

-- [Tensor] Tensor counterpart of Zeroes: an all-zero tensor of the given shape.
-- Unlike Eye_tensor, dims need not be a hypercube here -- an all-zero
-- tensor is well-defined for any shape.
local function Zeroes_tensor(dims)
    validate_tensor_dims(dims, "Zeroes_tensor", "Zeroes")
    return build_filled_tensor(dims, 1, 0)
end

-- [Matrix] If `seed` is given, reseeds Lua's global PRNG right before filling the
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

-- [Tensor]
local function build_random_tensor(dims, dim_idx, valor)
    local extent = dims[dim_idx]
    local t = {}
    if dim_idx == #dims then
        for i = 1, extent do
            t[i] = math.random(1, valor)
        end
        return t
    end
    for i = 1, extent do
        t[i] = build_random_tensor(dims, dim_idx + 1, valor)
    end
    return t
end

-- [Tensor] Tensor counterpart of Random_mat: same entry distribution (integers drawn
-- uniformly from math.random(1, valor)) and the same optional `seed` for
-- reproducibility, generalized to any shape with 3+ dimensions.
local function Random_tensor(dims, valor, seed)
    validate_tensor_dims(dims, "Random_tensor", "Random_mat")
    assert(type(valor) == "number" and valor >= 1, "Random_tensor: valor must be >= 1.")
    if seed ~= nil then
        assert(type(seed) == "number", "Random_tensor: seed must be a number.")
        math.randomseed(seed)
    end
    return build_random_tensor(dims, 1, valor)
end

-- [Matrix]
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

-- [Vector]
local function VTrace(matriz)
    local rows, cols = #matriz, #matriz[1]
    assert(rows == cols, "VTrace: matrix must be square.")
    local sum = 0
    for i = 1, rows do
        sum = sum + matriz[i][i]
    end
    return sum
end

-- [Vector] Vectors are plain flat tables {v1, v2, ...}, unlike matrices, which are
-- tables of row-tables. This is the bilinear (non-conjugated) dot product;
-- for the Hermitian inner product on complex vectors, conjugate one side
-- first with Conjugate().
local function Vdot(v1, v2)
    assert(#v1 == #v2, "Vdot: vectors must have the same length.")
    local sum = 0
    for i = 1, #v1 do
        sum = sum + v1[i] * v2[i]
    end
    return sum
end

-- [Vector]
local function VNorm(v, p)
    assert(type(v) == "table" and #v > 0, "VNorm: v must be a non-empty vector.")
    p = p or 2
    assert(p ~= 0, "VNorm: p must not be 0 (there is no p=0 'norm' expressible via sum(|vi|^p)^(1/p)).")
    local sum = 0
    for i = 1, #v do
        sum = sum + cabs(v[i]) ^ p
    end
    return sum ^ (1 / p)
end

-- [Vector] Returns a new vector scaled to unit length (VNorm(v, p) == 1), leaving `v`
-- itself untouched. Uses the same p as VNorm (default 2, Euclidean).
local function VNormalize(v, p)
    local norm = VNorm(v, p)
    assert(norm > EPSILON, "VNormalize: cannot normalize a zero vector.")
    local C = {}
    for i = 1, #v do
        C[i] = v[i] / norm
    end
    return C
end

-- [Determinant]
local function MDet(matrix)
    local rows, cols = #matrix, #matrix[1]
    assert(rows == cols, "MDet: matrix must be square.")

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

    if cabs(lu[rows][rows]) <= EPSILON then
        return 0
    end

    local det = sign
    for i = 1, rows do
        det = det * lu[i][i]
    end

    return det
end

-- ===========================================================================
-- Big_MDet: determinant of huge matrices via blocked LU decomposition with
-- partial pivoting, where the dominant O(n^3) step -- updating the trailing
-- submatrix after each panel -- is computed with Strassen's algorithm
-- (~O(n^2.807)) instead of naive triple-loop multiplication. This is the
-- standard "GEPP + Strassen" hybrid design real numerical libraries use: full
-- panels are still factored with ordinary partial pivoting (so it's exactly
-- as numerically stable as MDet -- there is no bad-pivot fragility here),
-- Strassen only accelerates the large rank-`block_size` update, which is
-- where virtually all the arithmetic actually happens for big n.
--
-- Below `block_size`, and for the recursion's own internal base case, this
-- just falls back to MDet's plain pivoted elimination / Mat_mul -- Strassen's
-- overhead (recursive submatrix allocation, 7 sub-multiplies + ~18 add/sub
-- passes per level) only pays off once n is large enough to amortize it.
-- Caveat worth being upfront about: in pure interpreted Lua, table-of-tables
-- allocation at every recursion level carries real constant-factor cost, so
-- the practical crossover where this beats MDet outright is workload- and
-- machine-dependent -- likely somewhere in the hundreds-to-low-thousands of
-- rows, not the dozens where Strassen starts winning in a compiled, cache-
-- friendly implementation. Tune `block_size` (default 64) if needed.
-- ===========================================================================

-- [Core] Shared by Big_MDet (Determinant) and FFT -- a genuine
-- cross-category dependency, unlike strassen_mul/strassen_mul_rect right
-- below, which only Big_MDet itself uses.
local function next_pow2(x)
    local m = 1
    while m < x do m = m * 2 end
    return m
end

-- [Determinant] Square Strassen multiply; n must already be a power of 2.
local function strassen_mul(A, B, n, threshold)
    if n <= threshold then
        return Mat_mul(A, B)
    end
    local h = n // 2
    local A11, A12, A21, A22 = {}, {}, {}, {}
    local B11, B12, B21, B22 = {}, {}, {}, {}
    for i = 1, h do
        A11[i], A12[i], A21[i], A22[i] = {}, {}, {}, {}
        B11[i], B12[i], B21[i], B22[i] = {}, {}, {}, {}
        for j = 1, h do
            A11[i][j] = A[i][j]
            A12[i][j] = A[i][j + h]
            A21[i][j] = A[i + h][j]
            A22[i][j] = A[i + h][j + h]
            B11[i][j] = B[i][j]
            B12[i][j] = B[i][j + h]
            B21[i][j] = B[i + h][j]
            B22[i][j] = B[i + h][j + h]
        end
    end

    local M1 = strassen_mul(Mat_sum(A11, A22), Mat_sum(B11, B22), h, threshold)
    local M2 = strassen_mul(Mat_sum(A21, A22), B11, h, threshold)
    local M3 = strassen_mul(A11, Mat_sub(B12, B22), h, threshold)
    local M4 = strassen_mul(A22, Mat_sub(B21, B11), h, threshold)
    local M5 = strassen_mul(Mat_sum(A11, A12), B22, h, threshold)
    local M6 = strassen_mul(Mat_sub(A21, A11), Mat_sum(B11, B12), h, threshold)
    local M7 = strassen_mul(Mat_sub(A12, A22), Mat_sum(B21, B22), h, threshold)

    local C11 = Mat_sum(Mat_sub(Mat_sum(M1, M4), M5), M7)
    local C12 = Mat_sum(M3, M5)
    local C21 = Mat_sum(M2, M4)
    local C22 = Mat_sum(Mat_sub(Mat_sum(M1, M3), M2), M6)

    local C = {}
    for i = 1, h do
        local row_top, row_bot = {}, {}
        for j = 1, h do
            row_top[j], row_top[j + h] = C11[i][j], C12[i][j]
            row_bot[j], row_bot[j + h] = C21[i][j], C22[i][j]
        end
        C[i], C[i + h] = row_top, row_bot
    end
    return C
end

-- [Determinant] Rectangular multiply (m x k) * (k x n) via square Strassen, using
-- zero-padding up to the next power of 2 and extracting the real m x n
-- corner of the result -- the zero-padded rows/columns only ever contribute
-- zero to the extracted region, so this is exact, not approximate.
local function strassen_mul_rect(A, m, k, B, k2, n, threshold)
    assert(k == k2, "strassen_mul_rect: inner dimensions must match.")
    local size = next_pow2(math.max(m, k, n))
    local Ap, Bp = {}, {}
    for i = 1, size do
        local row = {}
        for j = 1, size do
            row[j] = (i <= m and j <= k) and A[i][j] or 0
        end
        Ap[i] = row
    end
    for i = 1, size do
        local row = {}
        for j = 1, size do
            row[j] = (i <= k and j <= n) and B[i][j] or 0
        end
        Bp[i] = row
    end
    local Cp = strassen_mul(Ap, Bp, size, threshold)
    local C = {}
    for i = 1, m do
        local row = {}
        for j = 1, n do
            row[j] = Cp[i][j]
        end
        C[i] = row
    end
    return C
end

-- [Determinant]
local function Big_MDet(matrix, block_size)
    local n = #matrix
    assert(n == #matrix[1], "Big_MDet: matrix must be square.")
    block_size = block_size or 64
    -- A block_size <= 0 would make panel_end < col every iteration, so `col`
    -- never advances -- an infinite loop, not just a bad result.
    assert(type(block_size) == "number" and block_size >= 1 and block_size == math.floor(block_size),
        "Big_MDet: block_size must be a positive integer.")

    if n <= block_size then
        return MDet(matrix)
    end

    local lu = copy_matrix(matrix)
    local sign = 1
    local col = 1

    while col <= n do
        local panel_end = math.min(col + block_size - 1, n)

        -- Panel factorization: ordinary partial-pivoting elimination,
        -- touching only this panel's own columns (col..panel_end) for every
        -- row, panel or trailing alike -- deliberately NOT the full row
        -- width. This has to stay column-restricted here because a pivot
        -- swap on a *later* k within this same panel can still move a
        -- "trailing" row into panel territory; updating beyond-panel columns
        -- eagerly, keyed off a row's position at the time, silently uses the
        -- wrong (pre-swap) data for whichever row a later swap relocates.
        -- Row swaps themselves are still full-row (so A12 ends up correctly
        -- permuted), just not the column *updates*.
        for k = col, panel_end do
            local pivot_row = find_pivot_row(lu, k, n, k)
            if not pivot_row then
                return 0
            end
            if pivot_row ~= k then
                lu[k], lu[pivot_row] = lu[pivot_row], lu[k]
                sign = -sign
            end
            for i = k + 1, n do
                local factor = lu[i][k] / lu[k][k]
                for j = k + 1, panel_end do
                    lu[i][j] = lu[i][j] - factor * lu[k][j]
                end
                lu[i][k] = factor
            end
        end

        -- U12 = L11^-1 * A12, computed now that every row's final position
        -- within the panel is settled (no more swaps will touch rows
        -- col..panel_end after this point). A12 -- panel rows, columns
        -- beyond panel_end -- was left untouched above, so it still holds
        -- the original (but already correctly row-permuted) values; replay
        -- the panel's already-stored multipliers against just those columns,
        -- restricted to the `width`-many panel rows -- cheap regardless of n.
        for k = col, panel_end - 1 do
            for i = k + 1, panel_end do
                local factor = lu[i][k]
                for j = panel_end + 1, n do
                    lu[i][j] = lu[i][j] - factor * lu[k][j]
                end
            end
        end

        -- Trailing-submatrix update: instead of continuing the panel's
        -- row-by-row elimination out to column n (the O(n^3)-dominant part
        -- of plain Gaussian elimination), do it as one big
        -- Strassen-accelerated matrix multiply.
        if panel_end < n then
            local rest = n - panel_end
            local width = panel_end - col + 1
            local L21 = {}
            for i = 1, rest do
                local row = {}
                for j = 1, width do
                    row[j] = lu[panel_end + i][col + j - 1]
                end
                L21[i] = row
            end
            local U12 = {}
            for i = 1, width do
                local row = {}
                for j = 1, rest do
                    row[j] = lu[col + i - 1][panel_end + j]
                end
                U12[i] = row
            end
            local update = strassen_mul_rect(L21, rest, width, U12, width, rest, block_size)
            for i = 1, rest do
                local trailing_row = lu[panel_end + i]
                local update_row = update[i]
                for j = 1, rest do
                    trailing_row[panel_end + j] = trailing_row[panel_end + j] - update_row[j]
                end
            end
        end

        col = panel_end + 1
    end

    if cabs(lu[n][n]) <= EPSILON then
        return 0
    end

    local det = sign
    for i = 1, n do
        det = det * lu[i][i]
    end

    return det
end

-- ===========================================================================
-- FFT: radix-2 Cooley-Tukey Fast Fourier Transform, O(N log N) instead of
-- the O(N^2) a direct DFT sum would cost. Recursive divide-and-conquer,
-- same shape as strassen_mul above: split v into its even- and odd-indexed
-- samples, recursively transform each half, then combine with a single
-- "butterfly" pass. Every combine step multiplies by a twiddle factor built
-- as a genuine Complex value (via math.cos/math.sin), so -- exactly like
-- the rest of this module -- a real-valued input vector picks up Complex
-- entries "for free" through +/-/* the moment a twiddle factor touches it;
-- FFT never special-cases real vs. Complex input itself.
--
-- Only handles lengths that are an exact power of 2, reusing next_pow2 (see
-- Big_MDet above) to check that: the recursive halving here has to land on
-- exactly 1 at the bottom, which only happens for a power-of-2 start. This
-- is a real, honest limitation (not implementing Bluestein's algorithm or
-- similar to support arbitrary lengths) -- zero-pad v to the next power of
-- 2 yourself if you need to transform a different length; note that this
-- changes what is actually being transformed (a padded signal's spectrum is
-- not simply "the original spectrum at more points"), so it is left as an
-- explicit choice for the caller to make, not something FFT does silently.
-- ===========================================================================

-- [FFT] sign = -1 for the forward transform's e^{-i*2*pi*k*n/N} twiddle factors,
-- +1 for the inverse transform's e^{+i*2*pi*k*n/N} ones; IFFT additionally
-- divides every output by N afterwards, which this helper deliberately
-- does not do, since it would be wrong to apply twice across the recursion.
local function fft_recursive(v, sign)
    local n = #v
    if n == 1 then
        return { v[1] }
    end
    local half = n // 2
    local even, odd = {}, {}
    for i = 1, half do
        even[i] = v[2 * i - 1]
        odd[i] = v[2 * i]
    end
    local Feven = fft_recursive(even, sign)
    local Fodd = fft_recursive(odd, sign)
    local result = {}
    for k = 1, half do
        local angle = sign * 2 * math.pi * (k - 1) / n
        local w = Complex(math.cos(angle), math.sin(angle))
        local term = w * Fodd[k]
        result[k] = Feven[k] + term
        result[k + half] = Feven[k] - term
    end
    return result
end

-- [FFT] Forward FFT: X_k = sum_n v_n * e^(-i*2*pi*k*n/N). `v` may hold plain
-- numbers, Complex values, or a mix of both.
local function FFT(v)
    assert(type(v) == "table" and #v > 0, "FFT: v must be a non-empty vector.")
    local n = #v
    assert(next_pow2(n) == n,
        "FFT: v's length (" .. n .. ") must be a power of 2 -- the radix-2 Cooley-Tukey " ..
        "algorithm this implements only ever splits its input exactly in half; zero-pad v up " ..
        "to the next power of 2 first if it isn't already one.")
    return fft_recursive(v, -1)
end

-- [FFT] Inverse FFT: v_n = (1/N) * sum_k V_k * e^(+i*2*pi*k*n/N), recovering the
-- original (generally Complex-valued) sequence a forward FFT(v) produced.
-- Same power-of-2 length restriction as FFT, for the same reason.
local function IFFT(V)
    assert(type(V) == "table" and #V > 0, "IFFT: V must be a non-empty vector.")
    local n = #V
    assert(next_pow2(n) == n,
        "IFFT: V's length (" .. n .. ") must be a power of 2, for the same reason as FFT.")
    local result = fft_recursive(V, 1)
    for i = 1, n do
        result[i] = result[i] / n
    end
    return result
end

-- [Systems] Row-reduces a copy of the matrix (Gaussian elimination with partial
-- pivoting) and counts the pivots found. Works for non-square matrices too.
local function Rank(matriz)
    assert(type(matriz) == "table" and #matriz > 0 and type(matriz[1]) == "table" and #matriz[1] > 0,
        "Rank: matriz must be a non-empty matrix.")
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

-- [Systems] Gauss-Jordan elimination on [matrix | identity] with partial pivoting.
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
                if cabs(factor) > EPSILON then
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

-- [Systems] Solves Ax = b via Gaussian elimination with partial pivoting and back
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
            if cabs(factor) > EPSILON then
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

-- [Eigen] Modified Gram-Schmidt QR decomposition of a REAL square matrix (Q
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

-- [Eigen]
local function assert_real_matrix(matriz, fn_name)
    for i = 1, #matriz do
        for j = 1, #matriz[i] do
            assert(not is_complex(matriz[i][j]),
                fn_name .. ": Complex-entry input matrices are not supported (the QR " ..
                "iteration only runs in real arithmetic); Complex eigenvalues can still " ..
                "occur in the *output* for a real matrix with a complex-conjugate pair.")
        end
    end
end

-- [Eigen] Unshifted QR algorithm: iterating A <- R*Q (from A's QR decomposition)
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
                eigenvalues[i] = Complex(tr / 2, sq / 2)
                eigenvalues[i + 1] = Complex(tr / 2, -sq / 2)
            end
            i = i + 2
        else
            eigenvalues[i] = A[i][i]
            i = i + 1
        end
    end
    return eigenvalues
end

-- [Eigen] Same iteration as Eigenvalues, additionally accumulating the orthogonal
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

-- [Eigen] Modified Gram-Schmidt: takes a matrix whose ROWS are a (possibly linearly
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
            local coeff = Vdot(Conjugate(q), v)
            for d = 1, dim do
                v[d] = v[d] - coeff * q[d]
            end
        end
        if VNorm(v) > tol then
            basis[#basis + 1] = VNormalize(v)
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

-- [Eigen] Rayleigh quotient iteration: refines a user-supplied seed vector into an
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
    assert(VNorm(v0) > EPSILON, "Rayleigh: seed vector must not be the zero vector.")
    max_iter = max_iter or 100
    assert(type(max_iter) == "number" and max_iter >= 1 and max_iter == math.floor(max_iter),
        "Rayleigh: max_iter must be a positive integer.")
    tol = tol or 1e-12

    local v = VNormalize(v0)
    local mu = Vdot(Conjugate(v), mat_vec(matriz, v))
    local converged = false

    for _ = 1, max_iter do
        local shifted = Mat_sub(matriz, Scalar_mul(Eye(n), mu))
        local ok, v_next = pcall(Solve, shifted, v)
        if not ok then
            -- (A - mu*I) is numerically singular: mu is already essentially
            -- an eigenvalue -- exactly what convergence looks like, so treat
            -- it as success rather than propagating Solve's error.
            converged = true
            break
        end
        v = VNormalize(v_next)
        local mu_next = Vdot(Conjugate(v), mat_vec(matriz, v))
        local delta = mu_next - mu
        mu = mu_next
        if cabs(delta) < tol then
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
    Complex = Complex, Is_complex = is_complex, Conjugate = Conjugate,
    Show_matrix = Show_matrix, Show_tensor = Show_tensor, Matrix = Matrix, Tensor = Tensor,
    Sequence = Sequence, Show_table = Show_table,
    Mat_mul = Mat_mul, Mat_sum = Mat_sum, Mat_sub = Mat_sub, Scalar_mul = Scalar_mul,
    Eye = Eye, Eye_tensor = Eye_tensor, Zeroes = Zeroes, Zeroes_tensor = Zeroes_tensor,
    Random_mat = Random_mat, Random_tensor = Random_tensor, T = T, VTrace = VTrace,
    Vdot = Vdot, VNorm = VNorm, VNormalize = VNormalize, MDet = MDet, Big_MDet = Big_MDet, Rank = Rank, Inverse = Inverse,
    Solve = Solve, Eigenvalues = Eigenvalues, Eigenvectors = Eigenvectors,
    GRAM_SCH = GRAM_SCH, Rayleigh = Rayleigh, FFT = FFT, IFFT = IFFT,
}
