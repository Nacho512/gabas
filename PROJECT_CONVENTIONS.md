# GABAS project conventions

This file is the shared, tool-agnostic source of truth for how this codebase is
built. It applies regardless of which AI coding tool (Claude, Codex, or any
other) is working on it. Tool-specific files (`CLAUDE.md`, `AGENTS.md`, ...)
should point here rather than duplicate this content.

**Scope: this document covers `INTERNAL/LUA-PART-PROJECT/` only** -- the Lua
math library (`GABAS_MAT` and its `GABAS-*` domain projects). It does NOT
apply to `INTERNAL/LATEX-PART-PROJECT/` (the LaTeX side of `gabas`, e.g.
`gabas.sty`). The two are different languages with different conventions;
nothing here about naming, module structure, or Lua-specific numerical
discipline should be assumed to carry over to the LaTeX package. If/when the
LaTeX side accumulates its own conventions, they belong in a separate file
scoped to `INTERNAL/LATEX-PART-PROJECT/`, not here.

## Naming convention ("Nacho_case")

Every public function is `Capitalized_Snake_Case`: an initial capital, then
lowercase words joined by underscores at natural word boundaries.

Examples: `Derivative_at_point`, `Bisection`, `Compile_expression`,
`Binomial_coefficient`, `Decimal_to_hexadecimal`.

Private/internal helpers (never exported) are plain `lower_snake_case`.

## Project and module layout

Each top-level math domain is its own sibling project directory under
`INTERNAL/LUA-PART-PROJECT/GABAS-MAT/modules/`, e.g. `GABAS-LINAL`,
`GABAS-CALC-ONE-VAR`, `GABAS-DISCRETE-MATH`, `GABAS-COMBINATORICS`. Many more
exist as empty placeholders (`return {}`) for future domains -- check before
assuming a domain doesn't exist yet.

Inside each project:

```
GABAS-XXX/
  GABAS_XXX.lua        -- thin public aggregator (note: underscores, unlike the hyphenated directory name)
  modules/
    core.lua            -- input-validation plumbing, NEVER exported publicly
    <domain>.lua         -- the actual functions (may start as one file, split into
                            several once there's enough content to justify it --
                            "monolith first, submodules once earned")
  tests/
    run.lua              -- test harness, lists every test file to run
    test_<name>.lua       -- one file per feature/topic
```

The aggregator pattern (`GABAS_XXX.lua`): iterate a `submodules` list, `require`
each one, flatten its public exports into one table, with a guard that asserts
no two submodules export the same name. `core.lua` (and any other pure-plumbing
submodule, e.g. `dual.lua`, `quadrature.lua`, `tanh_sinh.lua`) is deliberately
**excluded** from that list -- it is internal infrastructure, not public API.

## The `core.lua` contract

Every project's `modules/core.lua`:

- `require`s `GABAS-LINAL.modules.core` directly (the true shared foundation)
  and re-exports its generic validators (`is_finite_number`, `is_integer`,
  `assert_positive_integer`, `assert_nonneg_integer`, `assert_nonblank_string`,
  ...) rather than re-deriving them.
- Adds only what's genuinely domain-specific on top (e.g. `assert_valid_base`
  in `GABAS-DISCRETE-MATH`).

Every public function validates its own inputs at the very top, via
`Core.assert_*` calls. Any invalid input stops execution **immediately** via
`assert`, with a message naming exactly which argument and function failed and
why (e.g. `"Bisection: a and b must be distinct."`). This is "Design by
Contract" in the formal sense: explicit preconditions that halt on violation
rather than letting bad input propagate into a confusing failure somewhere
downstream.

## Dependency direction (no circular `require`s)

- A cross-project `require` reaches directly into the specific submodule
  needed (e.g. `require("GABAS-LINAL.modules.complex")`), never through
  another project's own public aggregator -- keeps the dependency minimal and
  explicit, and avoids pulling in an entire public surface (with its
  duplicate-export check) just to reach one function.
- The whole dependency graph must stay a DAG. Lua's `require` handles a
  genuine cycle fragilely (a module can receive a half-built table from its
  still-loading partner).
- When two domains seem to need each other, get concrete before assuming a
  cycle is unavoidable:
  1. Check whether the need is really bidirectional at the level of one
     specific function -- it usually resolves to one direction once you're
     concrete about which function needs which.
  2. If genuinely mutual, extract the shared primitive into a third, lower
     module both depend on (e.g. `core.lua`/`dual.lua`).
  3. Split at a finer submodule granularity so only the specific submodule
     that truly needs the reverse dependency carries it -- a whole project
     "depending on" another because one of its many submodules needs it is a
     modeling mistake, not a real constraint.

## Numerical rigor

- Never improvise a numerical algorithm's exact structure or constants when a
  validated, well-known reference exists (e.g. the QK21 Gauss-Kronrod
  node/weight tables, copied from QUADPACK, not re-derived; Brent-Dekker's
  root finder, transcribed variable-for-variable from the classical
  reference algorithm). Transcription errors here are easy to make and easy
  to miss.
- Simple, well-known "textbook" algorithms (factorial, GCD, Fibonacci, base
  conversion, bisection) don't need an external reference -- implement
  directly and verify thoroughly with tests against independently-known
  values.
- When deriving a nontrivial formula from scratch (e.g. Filon quadrature's
  per-panel moment-integral coefficients), verify the derivation with a
  computer-algebra system (sympy or equivalent) before trusting it -- not
  from memory alone.
- Watch for catastrophic cancellation in a closed-form formula that has a
  removable singularity at some parameter value (e.g. an expression that is
  `0/0` in the limit `theta -> 0`). Use a Taylor-series branch below a
  numerically-verified threshold rather than trusting the closed form near
  the dangerous point.
- Never let a computation silently produce a wrong answer. Lua 5.4 integers
  wrap around silently on overflow -- guard against this explicitly (e.g.
  `Factorial`, `Binomial_coefficient`, `Fibonacci`) before it happens, not
  after.
- A numerical algorithm must never report a false "SUCCESS"/converged
  status. When a genuine limit is hit (max iterations, non-convergence, a
  mathematically nonexistent result), report an honest failure status
  instead of a plausible-looking wrong number.

## Verification discipline

- Never trust a claim about code behavior without actually running it
  against the real `lua5.4` interpreter.
- Every new function gets tests in its project's own `tests/` directory, run
  through that project's `tests/run.lua`.
- Tests should include: independently-known reference values (computed
  separately, e.g. via `mpmath`/Python -- never values the implementation
  itself produced), edge cases, input-validation checks (invalid
  types/out-of-range values correctly rejected with a specific message), and
  where practical, differential testing (cross-checking one implementation
  against an independent one -- e.g. `Brent` vs. `Bisection` on the same
  bracket).
- Run the full test suite of every project touched, not just the one with
  new code, before considering work done.

## Comments

- Default to **no comments**. Only add one when the WHY is genuinely
  non-obvious: a hidden constraint, a subtle invariant, a numerical hazard, a
  deliberate scope decision. Never explain WHAT the code does when the code
  (well-named identifiers) already says so.
- An explanatory comment authored by an AI coding assistant is prefixed with
  that assistant's own name, e.g. `-- Claude:` for Claude, `Codex:` for
  Codex -- so the authorship of any given piece of commentary stays
  traceable.

## Git / workflow discipline

- **Always ask for explicit confirmation before merging any PR to `main`.**
  This is a standing, unwavering policy -- never assumed, even after having
  been granted once before in an earlier session.
- Commit only a coherent, tested, working piece of work.
- Adding a new submodule or function always comes with its test file added to
  the relevant `tests/run.lua` list in the same change.

## Scope discipline

- Build the smallest verified piece first, then extend incrementally. Don't
  build a large speculative surface area in one shot.
- Don't add speculative generality or features beyond what's actually asked.
- Full symbolic manipulation / a CAS (parser -> mutable AST -> algebraic
  simplification -> symbolic differentiation/integration) is a deliberately
  **deferred**, later phase of this overall project -- mature the numeric
  layer first.
- When evaluating a large candidate feature list (e.g. an LLM-generated
  "complete" architecture survey), treat it as a wishlist/vision document,
  not a build plan. Extract only what's a natural, small, incremental
  extension of what already exists. Anything that is itself a
  library-sized undertaking (a full CAS, ODE solvers, FFT/transforms,
  probability distributions, plotting, certified interval arithmetic, ...)
  is its own future project, never a subsection folded into an existing one.

## Architectural decisions already made

- `GABAS_CALC_ONE_VAR.lua` is real-variable calculus only. Differentiating
  with respect to a *complex* variable itself (Jacobian/Wirtinger modes) is
  explicitly out of scope here, planned as a separate future module.
- Numeral-base conversion, GCD/LCM, and primality live in
  `GABAS-DISCRETE-MATH` (general discrete math).
- Pascal's triangle, binomial coefficients, and Fibonacci live in
  `GABAS-COMBINATORICS` (specifically combinatorial/recurrence-based), split
  out from `GABAS-DISCRETE-MATH` on purpose -- see that project's own file
  headers for the reasoning.
- Symbolic/CAS-style manipulation is its own future, later phase -- not to be
  threaded into the numeric modules prematurely, no matter how tempting a
  specific case looks.
