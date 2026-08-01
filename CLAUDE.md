# CLAUDE.md

This repo is `gabas`: a Lua math library (`GABAS_MAT`, under
`INTERNAL/LUA-PART-PROJECT/`) eventually meant to be embedded in a LaTeX
package (`INTERNAL/LATEX-PART-PROJECT/`). This session may be picking up
work from an earlier conversation that no longer has any context loaded.

## Read this first

**`PROJECT_CONVENTIONS.md`** (repo root) is the authoritative source for how
the Lua side of this project is built: naming convention, module/aggregator
layout, the `core.lua` input-validation contract, dependency-direction rules,
numerical rigor practices, verification/testing discipline, and the
architectural decisions already made. Read it before writing or changing any
Lua code. It explicitly does NOT cover the LaTeX side of the project.

## Orienting yourself

Each math domain lives in its own sibling directory under
`INTERNAL/LUA-PART-PROJECT/GABAS-MAT/modules/`, e.g. `GABAS-LINAL`,
`GABAS-CALC-ONE-VAR`, `GABAS-DISCRETE-MATH`, `GABAS-COMBINATORICS`. Many more
exist there as empty placeholders (`return {}`) for future domains -- check a
directory's actual contents before assuming it's unbuilt or already built.

To run a project's test suite (from anywhere, using an absolute path --
`tests/run.lua` locates its own root via `debug.getinfo`, so a relative `cd`
first is not required):

```
lua5.4 /home/user/gabas/INTERNAL/LUA-PART-PROJECT/GABAS-MAT/modules/<GABAS-PROJECT>/tests/run.lua
```

Run the full suite of every project touched by a change, not just the one
with new code, before considering the work done.

## Working style specific to this project

- Comments this session authors are prefixed `-- Claude:` (see
  `PROJECT_CONVENTIONS.md`'s comment-authorship convention -- each AI
  tool that has touched this repo gets its own prefix, e.g. `Codex:`).
- Always ask for explicit confirmation before merging any PR to `main` --
  standing policy, never assumed even if granted once before.
- The user (Nacho) is a lawyer, self-taught in math/programming, and
  explicitly wants real engineering rigor, not flattery or hand-waving --
  verify claims against the real `lua5.4` interpreter rather than asserting
  them, and say so plainly when something is a genuine limitation rather
  than smoothing it over.
