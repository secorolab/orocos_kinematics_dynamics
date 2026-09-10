# Plan 004: Retire the forked Vereshchagin variants and migrate motion-spec to the weighted solver

> **Executor instructions**: Follow this plan step by step. Run every
> verification command and confirm the expected result before moving to the
> next step. If anything in the "STOP conditions" section occurs, stop and
> report — do not improvise. When done, update the status row for this plan
> in `plans/README.md`.
>
> **Drift check (run first)**:
> `git diff --stat 4cadf46..HEAD -- orocos_kdl/src/`
> If any in-scope file changed since this plan was written, compare the
> "Current state" excerpts against the live code before proceeding; on a
> mismatch, treat it as a STOP condition.

## Status

- **Priority**: P2
- **Effort**: M
- **Risk**: MED
- **Depends on**: plans/001, plans/002, plans/003 — **all three must be
  merged into this fork's solver before starting**
- **Category**: tech-debt
- **Planned at**: commit `4cadf46`, 2026-08-07
- **Round**: 2 — **DEFERRED**. Deleting the forked classes breaks
  `domain_solver.stg`, which still names them, so repo A and repo B must merge
  in lockstep. Land plan 001 in the same round.

## Why this matters

This fork carries four solver files that upstream does not have, totalling
~2100 lines:

```
orocos_kdl/src/chainhdsolver_vereshchagin_fext.{hpp,cpp}
orocos_kdl/src/chainhdsolver_vereshchagin_fext_fixed_joint.{hpp,cpp}
orocos_kdl/src/chainhdsolver_vereshchagin_fixed_joint.{hpp,cpp}
```

They are copy-paste derivatives of `chainhdsolver_vereshchagin.cpp` with small
edits, and every upstream fix to the real solver has to be hand-ported into all
of them or silently miss them. Plans 002 and 003 put both capabilities —
fixed-joint support and driver prioritisation — into the one real solver, which
makes all four files deletable.

They are also **wrong** in a way that matters on hardware. The consuming
codegen runs the `_fext` solver as a separate pass and adds its torque to the
constrained pass (`motion-spec/src/motion_spec/templates/domain_solver.stg:551`).
That pass is supposed to carry the external-wrench contribution alone, but it
carries two more things:

- **Gravity.** Both solver instances are constructed with the same
  `root_acc = (0, 0, 9.81)`. Measured on a 3-link arm with **zero** external
  wrench, the `_fext` pass returns `(10.46, -20.14, 11.84)` Nm of pure gravity
  torque, which is then added on top of a constrained pass that already
  gravity-compensates.
- **Velocity-product bias.** `_fext` drops the `v x (H v)` term from `s.U` and
  the `PC` term from `totalBias`, but still carries `child.PC` into
  `s.R_tilde` (`chainhdsolver_vereshchagin_fext_fixed_joint.cpp:241`). With
  zero gravity and zero wrench but non-zero joint velocity it still emits
  `(0.67, -1.62, 1.30)` Nm.

Both numbers come from `plans/reference/achd_probe.cpp`, which you can re-run.

After this plan there is one solver, the arithmetic is a clean decomposition,
and the fork's diff against upstream is empty for `orocos_kdl/src/`.

## Current state

**This plan spans two repositories.** Paths below are relative to
`/home/batsy/work/ms/src/`.

### Repo A — `orocos_kinematics_dynamics` (this repo)

Fork-local changes against `origin/master`, from
`git diff --stat origin/master -- orocos_kdl/src/`:

```
 orocos_kdl/src/chainhdsolver_vereshchagin.cpp      |   4 +-
 orocos_kdl/src/chainhdsolver_vereshchagin.hpp      |   1 +
 orocos_kdl/src/chainhdsolver_vereshchagin_fext.cpp | 523 ++++++
 orocos_kdl/src/chainhdsolver_vereshchagin_fext.hpp | 535 ++++++
 ...chainhdsolver_vereshchagin_fext_fixed_joint.cpp | 518 ++++++
 ...chainhdsolver_vereshchagin_fext_fixed_joint.hpp | 535 ++++++
 .../src/chainhdsolver_vereshchagin_fixed_joint.cpp | 607 ++++++
 .../src/chainhdsolver_vereshchagin_fixed_joint.hpp | 536 ++++++
 orocos_kdl/src/joint.hpp                           |  10 +
 orocos_kdl/src/segment.hpp                         |   9 +
```

The 4-line change to `chainhdsolver_vereshchagin.cpp` is two separate things,
both from commit `42a9af9`:

```diff
+            E_input = s.E_tilde;
...
-        if (Sm(i) < 1e-14)
+        if (Sm(i) < 1e-8)
...
+    nu_sum += E_input.transpose() * acc;
```

The `E_input` term redefines `beta` to mean the *true* base-frame acceleration
instead of the gravity-offset one. It was measured to be **exactly** equivalent
to the caller passing `beta + alpha^T * root_acc` — max deviation `1.4e-14` in
`tau` and `q_dotdot` across a 12-pose sweep with non-zero velocity, feed-forward
torque and external wrench (`plans/reference/achd_beta_semantics.cpp`). It is a
change of input convention, not a fix, and it desynchronises `beta` from
`getTransformedLinkAcceleration()`. Plan 001 documents the convention upstream;
this plan moves the shift to the call site.

The SVD threshold change (`1e-14` → `1e-8`) is a separate, undocumented
behaviour change: it declares constraint directions singular a million times
more eagerly. **By explicit decision it is left as-is.** Do not revert it, do
not touch it in this plan. It is the one intentional remaining difference from
upstream in this file.

`joint.hpp` adds `Joint::setInertia(const double&)`; `segment.hpp` adds
`Segment::getMutableJoint()`. A search across `src/` found **no non-generated
consumer** of either. (`Segment::setInertia(const RigidBodyInertia&)` is a
different, pre-existing upstream method and is used by
`orocos_kdl/tests/solvertest.cpp:99` — do not confuse them.)

PyKDL binds only the fixed-joint variant, at
`python_orocos_kdl/PyKDL/kinfam.cpp:539-570`, with a Python test at
`python_orocos_kdl/tests/kinfamtest.py:158` registered at line 415.

### Repo B — `motion-spec`

`motion-spec/src/motion_spec/templates/domain_solver.stg` is the codegen. The
relevant fragments:

Construction of both solvers, lines 150–154:

```
<if(motion_driver.has_cartesian_force)>
        state.<solver.id>.f_ext = KDL::Wrenches(state.<solver.id>.num_segments);
        state.<solver.id>.achd_fext = std::make_unique\<KDL::ChainHdSolver_Vereshchagin_Fext_FixedJoint\>(*robot.<solver.id>.chain, state.<solver.id>.root_acc, state.<solver.id>.num_spatial_directions);
<endif>
        state.<solver.id>.achd_acc = std::make_unique\<KDL::ChainHdSolver_Vereshchagin_Fixed_Joint\>(*robot.<solver.id>.chain, state.<solver.id>.root_acc, state.<solver.id>.num_spatial_directions);
```

The extra pass, lines 500–518 (`solver-run-fext-pre-ACHD-robif2b`), which runs
the `_fext` solver with zeroed `alpha`, `beta` and `ff_torques`.

The wrench argument selector, lines 526–532 — note the `robif2b` backend passes
`f_ext_zero` to the constrained solver, so today the wrench reaches the
constrained solve on the `mj_kdl` backend only:

```
solver-achd-f_ext-arg-mj_kdl(solver, motion_driver) ::= <<
<if(motion_driver.has_cartesian_force)>state.<solver.id>.f_ext<else>f_ext_zero_<solver.id><endif>
>>

solver-achd-f_ext-arg-robif2b(solver, motion_driver) ::= <<
f_ext_zero_<solver.id>
>>
```

The torque sum, lines 549–555:

```
solver-run-finalize-robif2b(solver, motion_driver) ::= <<
<if(motion_driver.has_cartesian_force)>
    KDL::Add(tau_cstr_fext_<solver.id>, tau_ctrl_acc_<solver.id>, state.<solver.id>.tau_ctrl);
<else>
    state.<solver.id>.tau_ctrl = tau_ctrl_acc_<solver.id>;
<endif>
>>
```

Repo B convention, from `CLAUDE.md`: **all target-language text is rendered by
the template** — never build C++ strings in Python. Keep every code change
inside the `.stg`.

### Wrench provenance — this determines the migration per motion

The consuming models use `f_ext` in two different senses, and they migrate
differently:

- **Virtual wrench** (synthesised by a controller; does not physically exist,
  so its joint torque *must* be commanded). Example:
  `state.arm_solver_pick.f_ext[...] += shared.wrench_force_ctrl_pk_support_z`
  in the `pick_place_single` models.
- **Measured wrench** (read from an F/T sensor; physically exists, so
  commanding it would double it). Example: the `wrist_ft` sensor read in the
  `admittance_arc_single` models.

## Commands you will need

| Purpose | Command | Expected on success |
|---|---|---|
| Build KDL | `cd orocos_kdl && mkdir -p build && cd build && cmake -DENABLE_TESTS:BOOL=ON -DCMAKE_BUILD_TYPE=Release ./.. && make -j$(nproc)` | exit 0 |
| Test KDL | `cd orocos_kdl/build && make check` | exit 0 |
| Fork diff | `cd <repo A> && git diff --stat origin/master -- orocos_kdl/src/` | (see done criteria) |
| Build workspace pkg | `cbps <pkg>` from the workspace root | exit 0 |
| Re-run probes | see header comment in `plans/reference/achd_probe.cpp` | — |

`cbps` is this workspace's build wrapper; raw `colcon` fails under the agent
shell. Run examples in GUI mode, not `--headless`.

## Scope

**In scope**:
- Repo A: delete the six fork-local solver files; revert
  `chainhdsolver_vereshchagin.{hpp,cpp}` to upstream plus plans 001–003;
  `python_orocos_kdl/PyKDL/kinfam.cpp`;
  `python_orocos_kdl/tests/kinfamtest.py`
- Repo B: `motion-spec/src/motion_spec/templates/domain_solver.stg`

**Out of scope** (do NOT touch):
- `mj_kdl_wrapper/third_party/orocos_kinematics_dynamics/` — a **separate
  vendored copy** of KDL. It has its own lifecycle. Do not edit it here; note
  in your report whether it also carries the forked variants, because it will
  need the same migration later.
- `mj_kdl_wrapper/build_p013/` and any other `build*/` directory — build
  artifacts.
- Any file under `motion-spec-dsl/models/*/generated/` or
  `motion-spec-dsl/generation/` — these are **generated output**. They change
  by re-running codegen, never by hand.
- The `.motion` / DSL grammar and `ir_gen.py` — exposing weights as authorable
  DSL syntax is a deliberate follow-up, not this plan. This plan sets weights
  from the template with fixed values.

## Git workflow

- Repo A branch: `chore/retire-forked-vereshchagin-variants`
- Repo B branch: `feature/vereshchagin-driver-weights`
- Repo B is `secorolab/rec`-style: push to the existing PR branch, never create
  a side branch and never push to `main`.
- Do NOT push or open a PR unless the operator instructed it.

## Steps

### Step 1: Confirm the prerequisites actually landed

```
grep -n "setDriverWeights" orocos_kdl/src/chainhdsolver_vereshchagin.hpp
grep -c "Joint::Fixed" orocos_kdl/src/chainhdsolver_vereshchagin.cpp
```

**Verify**: the first returns a match; the second returns ≥ 6. If either
fails, plans 002/003 are not merged — STOP.

### Step 2: Add the caller-side gravity shift that replaces `E_input`

Plan 001 removes `E_input` from the solver. This step restores the equivalent
behaviour where it belongs — at the call site.

In repo B, `domain_solver.stg`, add the shift where `beta` is assembled: the
acceleration-energy setpoint passed to the solver becomes
`beta + alpha^T * root_acc`. Because `alpha` columns are unit vectors, this is
per-constraint-direction the dot product of that column with `root_acc`. Emit
it as a loop in the template next to where
`state.<solver.id>.acceleration_energy` is filled.

This is *exactly* equivalent to the removed `E_input` term — verified to
`1.4e-14` across a 12-pose sweep — so no model should change behaviour.

**Verify**: regenerate one model that uses a non-zero `root_acc`, and diff the
generated controller against the previous generation. The only difference
should be the added `beta` shift.

**Verify**: run that model in GUI; it must behave as before, not just compile.

### Step 3: Switch the codegen to one solver

In `domain_solver.stg`:

1. Replace both `ChainHdSolver_Vereshchagin_Fext_FixedJoint` and
   `ChainHdSolver_Vereshchagin_Fixed_Joint` with `ChainHdSolver_Vereshchagin`
   (lines 39, 152, 154). Drop the now-unused `achd_fext` member.
2. Delete the `solver-run-fext-pre-ACHD-robif2b` fragment (lines 500–518) and
   the `solver-run-fext-pre-impl` dispatch entry that selects it, so the
   `robif2b` backend no longer runs a second solve.
3. Make `solver-achd-f_ext-arg-robif2b` pass the real
   `state.<solver.id>.f_ext` instead of `f_ext_zero_<solver.id>`, so the
   wrench reaches the single constrained solve on both backends.
4. Replace `solver-run-finalize-robif2b`'s `KDL::Add(...)` with the plain
   assignment `state.<solver.id>.tau_ctrl = tau_ctrl_acc_<solver.id>;` for both
   branches.
5. Emit a `setDriverWeights` call in `solver-init-algorithm-ACHD`, after
   construction. Default both weight vectors to all-ones (classic behaviour).
   For a motion whose wrench is *measured* rather than virtual, emit
   `w_f_ext` all-zeros so the constraint does not fight the sensed force.
   **By explicit decision, hard-code all-ones for every model in this round.**
   The measured-vs-virtual wrench distinction is deferred; do not try to infer
   it per model. All-ones reproduces today's classic prioritisation, so no
   model changes behaviour from the weights.

**Verify**: `grep -rn "Fext_FixedJoint\|Vereshchagin_Fixed_Joint" motion-spec/src/`
→ no output.

### Step 4: Delete the forked solver files and their bindings

```
git rm orocos_kdl/src/chainhdsolver_vereshchagin_fext.hpp \
       orocos_kdl/src/chainhdsolver_vereshchagin_fext.cpp \
       orocos_kdl/src/chainhdsolver_vereshchagin_fext_fixed_joint.hpp \
       orocos_kdl/src/chainhdsolver_vereshchagin_fext_fixed_joint.cpp \
       orocos_kdl/src/chainhdsolver_vereshchagin_fixed_joint.hpp \
       orocos_kdl/src/chainhdsolver_vereshchagin_fixed_joint.cpp
```

`orocos_kdl/src/CMakeLists.txt` globs sources (`FILE( GLOB_RECURSE KDL_SRCS ...)`),
so no CMake edit is needed.

Then update PyKDL: in `python_orocos_kdl/PyKDL/kinfam.cpp:539-570` rebind the
block to `ChainHdSolver_Vereshchagin` (adding `setDriverWeights`), and rename
the Python test `testChainHdSolverVereshchaginFixedJoint` at
`python_orocos_kdl/tests/kinfamtest.py:158` and its registration at line 415.

**Verify**: `cd orocos_kdl/build && cmake ./.. && make -j$(nproc) && make check`
→ exit 0.

**Verify**: `grep -rn "Fixed_Joint\|_fext" orocos_kdl/src/ python_orocos_kdl/PyKDL/`
→ no output.

### Step 5: Evaluate the `joint.hpp` / `segment.hpp` additions

Search the whole workspace, excluding `build*/` and vendored third-party trees,
for `getMutableJoint` and `Joint::setInertia`. If genuinely unused, revert both
files to upstream. If used, keep them and note where — they would then need
their own upstream PR.

**Verify**: `git diff --stat origin/master -- orocos_kdl/src/joint.hpp orocos_kdl/src/segment.hpp`
→ empty, or a one-line justification recorded in your report.

### Step 6: Regenerate and run the affected models

Regenerate the `pick_place_single` and `admittance_arc_single` models and run
each in GUI mode to completion (`S_DONE`).

**Verify**: both reach `S_DONE`. Compare the logged joint torques against a
pre-migration run of the same model; the gravity contamination described in
"Why this matters" should disappear from the `robif2b` path.

## Test plan

- No new KDL unit tests — plans 002 and 003 supply those. This plan's tests are
  the existing suite plus the end-to-end model runs in Step 7.
- Re-run `plans/reference/achd_probe.cpp` after the migration: with the single
  weighted solver and `w_f_ext = 0`, the constrained solve must be independent
  of the wrench, and there must be no separate pass emitting gravity torque.
- Verification: `cd orocos_kdl/build && make check` → exit 0; both models reach
  `S_DONE`.

## Done criteria

ALL must hold:

- [ ] `git diff --stat origin/master -- orocos_kdl/src/` shows changes to
      `chainhdsolver_vereshchagin.{hpp,cpp}` **only** — the plans 002/003
      feature work plus the deliberately-kept `1e-8` SVD threshold — and no
      other file
- [ ] `cd orocos_kdl/build && make check` exits 0
- [ ] `grep -rn "Fext_FixedJoint\|Vereshchagin_Fixed_Joint\|E_input" orocos_kdl/src/ motion-spec/src/`
      returns no output
- [ ] `python tests/PyKDLtest.py` passes from `python_orocos_kdl/`
- [ ] `pick_place_single` and `admittance_arc_single` regenerate and reach `S_DONE` in GUI
- [ ] `plans/README.md` status row updated

## STOP conditions

Stop and report back (do not improvise) if:

- Plans 001/002/003 are not present in the solver (Step 1 fails).
- Removing the second solver pass changes a model's behaviour in a way you
  cannot attribute to the removed gravity/bias contamination. Report the
  measured torque difference per joint rather than re-tuning gains.
- You find a consumer of `getMutableJoint` or `Joint::setInertia` outside
  build artifacts. Keep those files and report.

Explicitly **not** STOP conditions — these are settled decisions, proceed:

- The IR does not distinguish measured from virtual wrenches. Expected.
  Hard-code all-ones weights (Step 3.5) and move on.
- The SVD threshold differs from upstream. Expected. Leave `1e-8` alone.

## Maintenance notes

- The real payoff is that `git diff origin/master -- orocos_kdl/src/` becomes
  reviewable. Keep it that way: the next capability belongs in the upstream
  solver behind a default-off switch, not in a fifth copied file.
- `mj_kdl_wrapper/third_party/orocos_kinematics_dynamics/` is a second vendored
  KDL copy and was deliberately left alone. Someone must check whether it also
  carries the forked variants and schedule the same migration.
- Deferred follow-ups: exposing the weights as authorable DSL syntax (grammar +
  `ir_gen.py` + template — the weights are per constraint direction, so they
  belong next to `alpha`/`beta` in the motion spec, not as solver constructor
  arguments); and binding `setDriverWeights` in PyKDL beyond the minimal
  rebinding in Step 5.
- A reviewer should scrutinise Step 3.3 — making the wrench reach the
  constrained solve on `robif2b` is the single behavioural change most likely
  to surprise on hardware. Bench-test before running on the arm.
