-- [Core]
-- Internal infrastructure shared by 2+ of the domain submodules (Matrix,
-- Tensor, Vector, Determinant, Systems, Eigen, FFT). Nothing in this file
-- is meant to be part of GABAS_LINAL.lua's public API -- every function
-- here is plumbing, not something an end user calls directly, which is why
-- GABAS_LINAL.lua's aggregator loop deliberately excludes "core" from the
-- list of submodules it flattens into the public table.
local Complex = require("GABAS-LINAL.modules.complex")

local MACHINE_EPSILON = 2.2204460492503131e-16
local DEFAULT_ABS_TOL = 1e-9
local DEFAULT_REL_TOL = 1e-9
-- Codex: backward-compatible internal name. New code should use the explicit
-- tolerance names or is_near_zero instead.
local EPSILON = DEFAULT_ABS_TOL

-- Magnitude that works on both plain numbers and Complex values -- used
-- everywhere a pivoting/convergence/near-zero check needs "how big is
-- this", regardless of whether the entries are real or complex.
local function cabs(x)
    if Complex.Is_complex(x) then
        return x:abs()
    end
    return math.abs(x)
end

-- Used by MDet/Big_MDet (Determinant), Rank (Systems), and
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

-- Searches rows [k, last_row] of column `col` for the entry with the
-- largest magnitude (classic partial pivoting: picking the largest
-- magnitude, not just the first non-zero one, keeps the elimination
-- numerically stable). Returns the row index, or nil if every candidate is
-- Codex: exact zero is singular here; scale-aware near-zero policy belongs
-- to the calling algorithm rather than an absolute pivot cutoff.
local function find_pivot_row(mat, k, last_row, col)
    local best_row, best_val = nil, 0
    for i = k, last_row do
        local v = cabs(mat[i][col])
        if v > best_val then
            best_row, best_val = i, v
        end
    end
    return best_row
end

-- Shared by Big_MDet (Determinant) and FFT -- a genuine cross-category
-- dependency, which is exactly why next_pow2 lives here instead of inside
-- either of those two submodules.
local function next_pow2(x)
    local m = 1
    while m < x do m = m * 2 end
    return m
end

-- True iff x is a plain Lua number that is neither NaN nor +-infinity.
-- (x == x is false exactly when x is NaN -- the standard portable check,
-- with no separate library call needed.)
local function is_finite_number(x)
    return type(x) == "number" and x == x and x ~= math.huge and x ~= -math.huge
end

-- True iff x is a finite real number OR a Complex value whose real and
-- imaginary parts are both finite -- "is this one legitimate scalar
-- entry", regardless of whether the surrounding matrix/vector is real or
-- complex. Every entry-validation loop in this module needs exactly this.
local function is_finite_scalar(x)
    if Complex.Is_complex(x) then
        return is_finite_number(x.re) and is_finite_number(x.im)
    end
    return is_finite_number(x)
end

-- True iff x is a finite number equal to its own floor -- "is this a
-- whole number", independent of sign (see is_positive_integer /
-- is_nonneg_integer below for the signed variants most call sites
-- actually want -- k, seeds, and iteration counts in this module are
-- never meant to be negative).
local function is_integer(x)
    return is_finite_number(x) and x == math.floor(x)
end

local function is_positive_number(x) return is_finite_number(x) and x > 0 end
local function is_nonneg_number(x) return is_finite_number(x) and x >= 0 end
local function is_positive_integer(x) return is_integer(x) and x >= 1 end
local function is_nonneg_integer(x) return is_integer(x) and x >= 0 end

-- True iff x is "numerically zero", real or complex, within `tol`
-- (defaulting to this module's general absolute tolerance) -- the boolean form of
-- the near-zero check find_pivot_row already does inline for pivoting.
local function is_zero(x, tol)
    return cabs(x) <= (tol or DEFAULT_ABS_TOL)
end

local function is_near_zero(x, scale, atol, rtol)
    scale = scale or 1
    atol = atol or DEFAULT_ABS_TOL
    rtol = rtol or DEFAULT_REL_TOL
    assert(is_finite_scalar(x), "Core.is_near_zero: x must be a finite scalar.")
    assert(is_nonneg_number(scale), "Core.is_near_zero: scale must be a finite nonnegative number.")
    assert(is_nonneg_number(atol), "Core.is_near_zero: atol must be a finite nonnegative number.")
    assert(is_nonneg_number(rtol), "Core.is_near_zero: rtol must be a finite nonnegative number.")
    return cabs(x) <= atol + rtol * scale
end

-- True iff x is a string with at least one non-whitespace character --
-- "not blank", for options like opts.method that must be a real,
-- meaningful name rather than "", " ", or nil.
local function is_nonblank_string(x)
    return type(x) == "string" and x:match("%S") ~= nil
end

-- assert_* wrap the predicates above with a standard message. `label`
-- should name both the calling function and the argument (e.g.
-- "SVD: opts.max_iter") so the error always says exactly what was wrong
-- and where -- the same discipline every hand-written assert elsewhere in
-- this module already follows, just no longer copy-pasted at every call
-- site that needs one of these checks.
local function assert_finite_number(x, label)
    assert(is_finite_number(x), label .. " must be a finite number.")
end

local function assert_finite_scalar(x, label)
    assert(is_finite_scalar(x), label .. " must be a finite number or a finite Complex value.")
end

local function assert_integer(x, label)
    assert(is_integer(x), label .. " must be an integer.")
end

local function assert_positive_number(x, label)
    assert(is_positive_number(x), label .. " must be a finite positive number.")
end

local function assert_nonneg_number(x, label)
    assert(is_nonneg_number(x), label .. " must be a finite nonnegative number.")
end

local function assert_positive_integer(x, label)
    assert(is_positive_integer(x), label .. " must be a positive integer.")
end

local function assert_nonneg_integer(x, label)
    assert(is_nonneg_integer(x), label .. " must be a nonnegative integer.")
end

local function assert_nonblank_string(x, label)
    assert(is_nonblank_string(x), label .. " must be a non-blank string.")
end

-- Generic enum check: x must equal one of the entries in `allowed` (a
-- plain array, e.g. {"QR", "AUTO"}). One assert instead of every call site
-- hand-rolling its own chain of `x == "A" or x == "B" or ...`.
local function assert_one_of(x, allowed, label)
    for i = 1, #allowed do
        if x == allowed[i] then return end
    end
    local names = {}
    for i = 1, #allowed do names[i] = tostring(allowed[i]) end
    assert(false, label .. " must be one of: " .. table.concat(names, ", ") .. ".")
end

-- kind (string) -> the assert_* above it dispatches to. "one_of" isn't
-- listed here -- it needs the extra `allowed` list argument, so validate()
-- below special-cases it directly instead of trying to fit it into this
-- same one-argument shape.
local VALIDATORS = {
    finite_number = assert_finite_number,
    finite_scalar = assert_finite_scalar,
    integer = assert_integer,
    positive_number = assert_positive_number,
    nonneg_number = assert_nonneg_number,
    positive_integer = assert_positive_integer,
    nonneg_integer = assert_nonneg_integer,
    nonblank_string = assert_nonblank_string,
}

-- Master validator: a single entry point over every assert_* above,
-- selected by name instead of having to remember which specific function
-- to call. `extra` is only used by kind = "one_of", where it's the
-- allowed-values array: validate(method, "one_of", "SVD_reduced: opts.method",
-- {"QR", "AUTO"}). Everything this dispatches to is still directly
-- callable on its own too (Core.assert_positive_integer, etc.) -- validate()
-- is a convenience layered on top, not a replacement for them.
local function validate(x, kind, label, extra)
    if kind == "one_of" then
        assert_one_of(x, extra, label)
        return
    end
    local fn = VALIDATORS[kind]
    assert(fn, "Core.validate: unknown validation kind \"" .. tostring(kind) .. "\".")
    fn(x, label)
end

local function assert_vector(v, label, opts)
    opts = opts or {}
    assert(type(opts) == "table", "Core.assert_vector: opts must be a table.")
    label = label or "vector"
    assert(type(v) == "table" and (opts.allow_empty or #v > 0),
        label .. " must be " .. (opts.allow_empty and "a vector." or "a non-empty vector."))
    if opts.length ~= nil then
        assert(#v == opts.length, label .. " must have length " .. opts.length .. ".")
    end
    for i = 1, #v do
        assert_finite_scalar(v[i], label .. "[" .. i .. "]")
        if opts.real then
            assert(not Complex.Is_complex(v[i]), label .. "[" .. i .. "] must be a finite real number.")
        end
    end
    return #v
end

local function assert_matrix(matrix, label, opts)
    opts = opts or {}
    assert(type(opts) == "table", "Core.assert_matrix: opts must be a table.")
    label = label or "matrix"
    assert(type(matrix) == "table" and #matrix > 0 and type(matrix[1]) == "table" and #matrix[1] > 0,
        label .. " must be a non-empty matrix (a table of non-empty row-tables).")
    local rows, cols = #matrix, #matrix[1]
    for i = 1, rows do
        assert(type(matrix[i]) == "table" and #matrix[i] == cols,
            label .. " must be rectangular -- every row must have length " .. cols .. ".")
        for j = 1, cols do
            assert_finite_scalar(matrix[i][j], label .. "[" .. i .. "][" .. j .. "]")
            if opts.real then
                assert(not Complex.Is_complex(matrix[i][j]),
                    label .. "[" .. i .. "][" .. j .. "] must be a finite real number.")
            end
        end
    end
    if opts.square then
        assert(rows == cols, label .. " must be square.")
    end
    return rows, cols
end

return {
    MACHINE_EPSILON = MACHINE_EPSILON,
    DEFAULT_ABS_TOL = DEFAULT_ABS_TOL,
    DEFAULT_REL_TOL = DEFAULT_REL_TOL,
    EPSILON = EPSILON,
    cabs = cabs,
    copy_matrix = copy_matrix,
    find_pivot_row = find_pivot_row,
    next_pow2 = next_pow2,
    is_finite_number = is_finite_number,
    is_finite_scalar = is_finite_scalar,
    is_integer = is_integer,
    is_positive_number = is_positive_number,
    is_nonneg_number = is_nonneg_number,
    is_positive_integer = is_positive_integer,
    is_nonneg_integer = is_nonneg_integer,
    is_zero = is_zero,
    is_near_zero = is_near_zero,
    is_nonblank_string = is_nonblank_string,
    assert_finite_number = assert_finite_number,
    assert_finite_scalar = assert_finite_scalar,
    assert_integer = assert_integer,
    assert_positive_number = assert_positive_number,
    assert_nonneg_number = assert_nonneg_number,
    assert_positive_integer = assert_positive_integer,
    assert_nonneg_integer = assert_nonneg_integer,
    assert_nonblank_string = assert_nonblank_string,
    assert_one_of = assert_one_of,
    validate = validate,
    assert_vector = assert_vector,
    assert_matrix = assert_matrix,
}
