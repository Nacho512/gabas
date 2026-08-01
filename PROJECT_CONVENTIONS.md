# GABAS project conventions

This file is the shared, tool-agnostic source of truth for how this codebase
is built. It applies regardless of which AI coding tool (Claude, Codex, or any
other) is working on it. Tool-specific files (`CLAUDE.md`, `AGENTS.md`, ...)
should point here rather than duplicate this content.

**Scope: this document covers the whole `gabas` project** -- both
`INTERNAL/LUA-PART-PROJECT/` (the Lua math library, `GABAS_MAT` and its
`GABAS-*` domain projects) and `INTERNAL/LATEX-PART-PROJECT/` (the LaTeX
package, e.g. `gabas.sty`). It's organized in three parts:

1. **General conventions** -- apply to the project as a whole, regardless of
   language or which side you're working on.
2. **Lua-specific conventions** -- apply only to `INTERNAL/LUA-PART-PROJECT/`.
3. **LaTeX-specific conventions** -- apply only to
   `INTERNAL/LATEX-PART-PROJECT/`.

The two sides are different languages with very different conventions --
notably, they disagree on naming case by design (see each section). Nothing
in the Lua-specific or LaTeX-specific parts should be assumed to carry over to
the other side; only Part 1 does.

---

## Part 1: General conventions

These apply everywhere in the project, on both the Lua and LaTeX sides.

### Comments

- Default to **no comments**. Only add one when the WHY is genuinely
  non-obvious: a hidden constraint, a subtle invariant, a numerical hazard, a
  deliberate scope decision. Never explain WHAT the code does when the code
  (well-named identifiers) already says so.
- An explanatory comment authored by an AI coding assistant is prefixed with
  that assistant's own name, e.g. `-- Claude:` for Claude, `Codex:` for
  Codex -- so the authorship of any given piece of commentary stays
  traceable.

### Verification discipline

- Never trust a claim about code behavior without actually running it, on
  whichever real toolchain applies to that side of the project (the
  `lua5.4` interpreter for the Lua side; actual engine compilation for the
  LaTeX side). Confident-sounding prose is not a substitute for output from
  the real tool.
- Say so plainly when something is a genuine limitation rather than smoothing
  it over.

### Git / workflow discipline

- **Always ask for explicit confirmation before merging any PR to `main`.**
  This is a standing, unwavering policy -- never assumed, even after having
  been granted once before in an earlier session.
- Commit only a coherent, tested, working piece of work.

### Scope discipline

- Build the smallest verified piece first, then extend incrementally. Don't
  build a large speculative surface area in one shot.
- Don't add speculative generality or features beyond what's actually asked.
- When evaluating a large candidate feature list (e.g. an LLM-generated
  "complete" architecture survey), treat it as a wishlist/vision document,
  not a build plan. Extract only what's a natural, small, incremental
  extension of what already exists. Anything that is itself a
  library-sized undertaking is its own future project, never a subsection
  folded into an existing one.

---

## Part 2: Lua-specific conventions

**Scope: `INTERNAL/LUA-PART-PROJECT/` only.**

### Naming convention ("Nacho_case")

Every public function is `Nacho_case_innstance`: an initial capital, then
lowercase words joined by underscores at natural word boundaries.

Examples: `Derivative_at_point`, `Bisection`, `Compile_expression`,
`Binomial_coefficient`, `Decimal_to_hexadecimal`.

Private/internal helpers (never exported) are plain `lower_snake_case`.

(Note: this is the opposite convention from the LaTeX side, which prefers
lowercase names throughout -- see Part 3. Do not mix the two.)

### Project and module layout

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

### The `core.lua` contract

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

### Dependency direction (no circular `require`s)

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

### Numerical rigor

- Never improvise a numerical algorithm's exact structure or constants when a
  validated, well-known reference exists (e.g. the QK21 Gauss-Kronrod
  node/weight tables, copied from QUADPACK, not re-derived; Brent-Dekker's
  root finder, transcribed variable-for-variable from the classical
  reference algorithm; the Lanczos coefficients for Gamma; Miller's
  backward-recurrence algorithm for Bessel_j). Transcription errors here are
  easy to make and easy to miss.
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
- Compute in LOG-SPACE FIRST whenever an intermediate value could
  over/underflow even though the final result is representable (e.g.
  `Log_gamma` before `Gamma`, `Log_beta` before `Beta`, the first term of
  `Modified_bessel_i`'s series via `Log_gamma(n+1)` instead of forming `n!`
  and `half_x^n` as separate numbers). This exact class of bug has been
  caught empirically, repeatedly, across unrelated functions -- always
  worth checking for on any new function whose inputs can get large.
- Never let a computation silently produce a wrong answer. Lua 5.4 integers
  wrap around silently on overflow -- guard against this explicitly (e.g.
  `Factorial`, `Binomial_coefficient`, `Fibonacci`, and any accumulator like
  a running factorial inside a series) before it happens, not after. When in
  doubt, seed an accumulator as a float literal (`1.0`, not `1`) so Lua
  arithmetic stays in float space throughout.
- A numerical algorithm must never report a false "SUCCESS"/converged
  status. When a genuine limit is hit (max iterations, non-convergence, a
  mathematically nonexistent result), report an honest failure status
  instead of a plausible-looking wrong number.
- When a numerical method has a verified-but-limited domain of accuracy (not
  a theoretical convergence limit but an empirically-determined precision
  boundary), enforce that domain explicitly via input validation rather than
  silently extrapolating past it -- unless/until a more general method
  (e.g. a different algorithm with no such boundary) replaces it.

### Verification discipline (Lua-specific mechanics)

Builds on the general verification principle in Part 1:

- Every new function gets tests in its project's own `tests/` directory, run
  through that project's `tests/run.lua`.
- Tests should include: independently-known reference values (computed
  separately, e.g. via `mpmath`/Python -- never values the implementation
  itself produced), edge cases, input-validation checks (invalid
  types/out-of-range values correctly rejected with a specific message), and
  where practical, differential testing (cross-checking one implementation
  against an independent one -- e.g. `Brent` vs. `Bisection` on the same
  bracket).
- A bug caught this way gets a dedicated regression test using the exact
  case that exposed it, not just a fix.
- Run the full test suite of every project touched, not just the one with
  new code, before considering work done.
- Adding a new submodule or function always comes with its test file added
  to the relevant `tests/run.lua` list in the same change.
- Full symbolic manipulation / a CAS (parser -> mutable AST -> algebraic
  simplification -> symbolic differentiation/integration) is a deliberately
  **deferred**, later phase of this overall project -- mature the numeric
  layer first.

### Architectural decisions already made

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

---

## Part 3: LaTeX-specific conventions

**Scope: `INTERNAL/LATEX-PART-PROJECT/` only.**

1. Write the LaTeX layer engine-independently whenever practical and
   possible; where engine-specific behavior is necessary, design primarily
   for LuaLaTeX.
2. Keep major capabilities modular. OCR, drawing recognition, Python
   integrations, Lua mathematics, and similar systems should live in
   dedicated modules loaded by `gabas.sty`.
3. Avoid forcing users to switch engines for different features.
   Engine-specific implementations should be selected internally through
   explicit capability detection.
4. Validate every public argument before filesystem access, PDF inspection,
   arithmetic, or rendering. Values must match their documented types
   exactly; invalid integers must never be silently rounded or truncated.
5. Distinguish recoverable problems from nonrecoverable ones. Recoverable
   conditions receive warnings and controlled fallbacks; invalid states that
   make execution unsafe or ambiguous receive fatal GABAS-specific errors.
6. Preserve diagnostics, compilation counters, manifests, routing
   information, and draft-mode recommendations. These mechanisms should be
   improved progressively rather than removed.
7. Keep obsolete or incorrect code as comments instead of deleting it,
   clearly marking replacement code and the reason for the correction.
8. Prefer lowercase names for new files, modules, internal identifiers, and
   generated artifacts whenever compatibility permits.
9. Avoid using names for functions, variables, etc. which contain uppercase
   letters as much as possible.
10. Follow the expl3 syntax-boundary rule below for every expl3
    implementation section.

(Note: rules 8-9 are the opposite convention from the Lua side, which
requires `Capitalized_Snake_Case` for every public function -- see Part 2.
Do not mix the two.)

### expl3 syntax-boundary rule

Every contiguous expl3 implementation section must begin with its own
explicit `\ExplSyntaxOn` and end with its own explicit `\ExplSyntaxOff`.

A contiguous expl3 implementation section is one uninterrupted region of
source code that collectively implements a single public command, internal
subsystem, or closely related family of commands using expl3 syntax.

Such a section may contain:

- expl3 variable declarations;
- message declarations;
- key definitions;
- private validation, parsing, filesystem, logging, calculation, or
  rendering helpers;
- hooks required exclusively by that subsystem;
- the public command definitions that invoke those helpers.

These elements belong in the same section only when they share one clear
responsibility and are maintained as one logical implementation. For
example, `\insertnacho`, its keys, messages, variables, validators,
PDF-routing helpers, logger, renderer, and public-command definition
constitute one contiguous implementation section.

The section must end before encountering:

- unrelated package functionality;
- a different public subsystem with independent responsibilities;
- substantial traditional LaTeX2e code that does not require expl3 syntax;
- document-facing declarations unrelated to the current subsystem;
- another implementation that could be removed without affecting the
  current subsystem.

#### Required structure

```latex
% Image-management subsystem
\ExplSyntaxOn

% Variables, messages, keys, helpers, hooks, and public commands
% belonging exclusively to this subsystem.

\ExplSyntaxOff
```

#### Additional requirements

1. `\ExplSyntaxOn` must never be left active merely for convenience.
2. An unrelated subsystem must receive a separate syntax-boundary pair.
3. Private helpers must remain beside the subsystem they serve whenever
   practical.
4. A section must not be split merely because it contains several helper
   functions.
5. Two implementations must not be combined merely because both use expl3.
6. Ordinary LaTeX code may occur inside the section when it is an integral
   part of that expl3 subsystem.
7. The delimiters provide syntax-state boundaries, not TeX grouping,
   variable locality, namespaces, or runtime isolation.
8. If responsibility is ambiguous, use the smaller coherent section and
   close `\ExplSyntaxOff` earlier.
