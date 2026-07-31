-- LATER ON I SHALL COMPLETE THE Á LA COBOL DOCUMENTATION
-- NACHI
-- 2026/07/30,10:10

-- ===========================================================================
--
-- GABAS_LINAL: thin aggregator over the modules/ submodules below. Splitting
-- the actual implementation into per-category files keeps each one small
-- enough to read and edit in isolation; this file's only job is to require
-- every domain submodule and flatten their public functions into one table,
-- so calling code is completely unaffected: `GL.MDet(A)`, `GL.FFT(v)`, etc.
-- work exactly as they did when everything lived in a single file.
--
--   Core                    -- modules/core.lua                      -- internal plumbing, not public
--   Domain                  -- modules/domain.lua                    -- implements scalar-domain abstractions, coercion, and domain capabilities
--   Complex                 -- modules/complex.lua                   -- provides complex numbers definition, manipulation and operations; necessary for many linear algebra operatios
--   Scalar                  -- modules/scalar.lua                    -- provides generic scalar operations routed through scalar domains
--   Matrix                  -- modules/matrix.lua                    -- implements matrix construction, operations, access, arithmetic, structural editing, and predicates
--   Tensor                  -- modules/tensor.lua                    —- implements tensors, tensor products, contractions, and exterior operations.
--   Vector                  -- modules/vector.lua                    -- implements vectors, vector arithmetic, products, geometry, and utilities
--   Permutation            -- modules/permutation.lua               -- implements permutation objects and their matrix actions.
--   Systems                 -- modules/systems.lua                   —- implements exact linear-system construction, classification, and complete solution sets
--   Eigen                   -- modules/eigen.lua                     -- Implements exact and formal eigenstructure.
--   NLA                     -- modules/nla.lua                       -- numerical linear algebra, a very useful and practical subject
--   SIMPLEX                 -- modules/simplex.lua                   -- esential for process engineers and managers
--   rowreduce           -- modules/rowreduce.lua                     — Implements exact row operations, REF, RREF, pivots, and proof traces
--   space               -- modules/space.lua                         — Implements vector spaces and global finite-dimensional space operations.
--   subspace            -- modules/subspace.lua                      — Implements spans, intersections, complements, quotients, and membership.
--   basis               -- modules/basis.lua                         — Implements basis validation, extension, reduction, coordinates, and transitions.
--   linear_map          -- modules/linear_map.lua                  — Implements linear transformations and their algebra.
--   determinant              -- modules/determinant.lua               — Implements determinants, minors, cofactors, adjugates, and determinant derivations.
--   polynomial             -- modules/polynomial.lua                 — Implements the polynomial objects required by characteristic and minimal-polynomial operations.
--   canonical             -- modules/canonical.lua                    — Implements canonical forms and associated transformations.
--   inner_product         -- modules/inner_product.lua                   — Implements inner-product spaces, orthogonality, projections, and adjoints.
--   dual                   -- modules/dual.lua                             — Implements dual spaces, functionals, annihilators, and transpose maps.
--   forms                    -- modules/forms.lua                          — Implements bilinear, sesquilinear, Hermitian, and quadratic forms.
--   affine                  -- modules/affine.lua                          — Implements affine spaces, affine combinations, and affine maps.
--   modules              -- modules/modules.lua                            — Implements the advanced module-over-rings subsystem.
--   symbolic                -- modules/symbolic.lua                        — Connects exact symbolic simplification, substitution, and assumptions.
--   parser                  -- modules/parser.lua                          — Implements controlled input grammars.
--   latex                 -- modules/latex.lua                             — Converts mathematical objects and proof traces into LaTeX.
--   serialize              -- modules/serialize.lua                        — Implements JSON, CSV, Matrix Market, and internal persistence.
--   diagnostics            -- modules/diagnostics.lua                      — Defines structured errors, warnings, validation messages, and result diagnostics.
--   testing               -- modules/testing.lua                           — Defines assertions, generators, and verification functions
--   FFT                   -- modules/fft.lua
--
-- ===========================================================================

local LA = {}

-- "core" is deliberately NOT in this list: modules/core.lua holds internal
-- plumbing (EPSILON, cabs, copy_matrix, find_pivot_row, next_pow2) that the
-- domain submodules below require directly from each other, not something
-- meant for the end user. If it were merged in here too, this same loop
-- would silently promote it to public API.

local submodules = {
    "complex",
    "matrix",
    "determinant",
    "tensor",
    "vector",
    "systems",
    "eigen",
    "fft",
    "nla",
    "simplex",
    "permutation",
    "scalar",
    "rowreduce",
    "space",
    "subspace",
    "basis",
    "linear_map",
    "polynomial",
    "canonical",
    "inner_product",
    "dual",
    "forms",
    "affine",
    "modules",
    "symbolic",
    "parser",
    "latex",
    "serialize",
    "diagnostics",
    "testing",
}

for _, name in ipairs(submodules) do
    local mod = require("GABAS-LINAL.modules." .. name)
    for key, value in pairs(mod) do
        assert(LA[key] == nil, "GABAS_LINAL: duplicate public export '" .. key .. "' from module " .. name .. ".")
        LA[key] = value
    end
end

return LA




 -- B"H.
