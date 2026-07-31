-- [Tensor]
-- n-dimensional generalization of the Matrix family. Requires Complex only
-- for Is_complex, used by tensor_depth to detect where a tensor's nesting
-- bottoms out at a leaf value.
local Complex = require("GABAS-LINAL.modules.complex")
local Core = require("GABAS-LINAL.modules.core")

-- NOTE (Nacho): FFT will likely become indispensable here once tensor-domain
-- operations are added, but nothing in this file needs it yet -- not
-- required until something in this module actually calls into it.

-- Shared guard for the whole *_tensor family (Tensor, Zeroes_tensor,
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
        Core.assert_positive_integer(extent, fn_name .. ": dims[" .. i .. "]")
    end
end

-- Builds a tensor of shape dims[dim_idx..#dims], every entry set to
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

-- Counts how many levels of table-nesting `t` has before hitting a leaf
-- (a plain number or a Complex value) -- i.e. how many dimensions a tensor
-- built by Tensor/Zeroes_tensor/Eye_tensor/Random_tensor actually has.
-- Only ever walks down through t[1], t[1][1], ...: like the rest of the
-- module, a ragged (non-rectangular) tensor is assumed not to occur.
local function tensor_depth(t)
    if type(t) ~= "table" or Complex.Is_Complex(t) then
        return 0
    end
    return 1 + tensor_depth(t[1])
end

-- Prints one (depth-2)-dimensional matrix "slice" of a tensor per call,
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

-- Tensor counterpart of Show_matrix: prints every 2-D slice of a tensor of
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

-- Recursively fills the innermost-to-outermost nesting of a tensor, one
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

-- Tensor generalizes Matrix from 2 dimensions to any number of them: `data`
-- stays a flat vector (exactly like Matrix's `data`), and `dims` is a table
-- listing EVERY dimension's size explicitly -- e.g. {2, 3, 4} for a 2x3x4
-- tensor -- rather than trying to infer all-but-one dimension from a single
-- count the way Matrix infers columns from num_filas. That inference does
-- not generalize past 2 dimensions (there is no unique way to split a
-- leftover element count across 2+ unknown dimensions), so for a tensor
-- every dimension size must be given up front. The result is nested tables
-- #dims levels deep: a 3-D tensor T is indexed T[i][j][k].
local function Tensor(data, dims)
    Core.assert_vector(data, "Tensor: data")
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

-- Builds levels dim_idx..#dims of the hyperdiagonal, given that every index
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

-- Tensor counterpart of Eye. A matrix's identity property (Eye(n) is the
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

-- Tensor counterpart of Zeroes: an all-zero tensor of the given shape.
-- Unlike Eye_tensor, dims need not be a hypercube here -- an all-zero
-- tensor is well-defined for any shape.
local function Zeroes_tensor(dims)
    validate_tensor_dims(dims, "Zeroes_tensor", "Zeroes")
    return build_filled_tensor(dims, 1, 0)
end

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

-- Tensor counterpart of Random_mat: same entry distribution (integers drawn
-- uniformly from math.random(1, valor)) and the same optional `seed` for
-- reproducibility, generalized to any shape with 3+ dimensions.
local function Random_tensor(dims, valor, seed)
    validate_tensor_dims(dims, "Random_tensor", "Random_mat")
    Core.assert_positive_integer(valor, "Random_tensor: valor")
    if seed ~= nil then
        Core.assert_integer(seed, "Random_tensor: seed")
        math.randomseed(seed)
    end
    return build_random_tensor(dims, 1, valor)
end

return {
    Tensor = Tensor,
    Zeroes_Tensor = Zeroes_tensor,
    Eye_Tensor = Eye_tensor,
    Random_Tensor = Random_tensor,
    Show_Tensor = Show_tensor,
}
