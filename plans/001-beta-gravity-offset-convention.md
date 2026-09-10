# Plan 001: Document and pin the `beta` gravity-offset convention in ChainHdSolver_Vereshchagin

> **Executor instructions**: Follow this plan step by step. Run every
> verification command and confirm the expected result before moving to the
> next step. If anything in the "STOP conditions" section occurs, stop and
> report — do not improvise. When done, update the status row for this plan
> in `plans/README.md`.
>
> **Drift check (run first)**:
> `git diff --stat 4cadf46..HEAD -- orocos_kdl/src/chainhdsolver_vereshchagin.hpp orocos_kdl/src/chainhdsolver_vereshchagin.cpp orocos_kdl/tests/solvertest.cpp`
> If any in-scope file changed since this plan was written, compare the
> "Current state" excerpts against the live code before proceeding; on a
> mismatch, treat it as a STOP condition.

## Status

- **Priority**: P1
- **Effort**: S
- **Risk**: LOW
- **Depends on**: none
- **Category**: docs
- **Planned at**: commit `4cadf46`, 2026-08-07
- **Round**: 2 — **DEFERRED**, must merge in lockstep with the `motion-spec`
  template PR

> **Why this is round 2 and not round 1.** Step 1 removes the fork's `E_input`
> term, which silently redefines what `beta` means for every generated
> controller. Shipping it without the matching caller-side shift in
> `motion-spec/src/motion_spec/templates/domain_solver.stg` would change the
> behaviour of every existing model with no compile error to warn anyone. Do
> not start this plan until the motion-spec change is ready to land with it.
>
> **Commit split.** Two commits, in this order:
> 1. `chore:` remove `E_input` — fork-only, **not** cherry-picked upstream.
> 2. `docs:` the convention documentation + test — upstream-clean, this is the
>    commit that becomes the upstream PR.
>
> The order is forced: the Step 4 test cannot pass while `E_input` is present.

## Why this matters

`ChainHdSolver_Vereshchagin` absorbs gravity by setting the root acceleration
to `-a_g` (Featherstone's trick, RBDA §5.3). A consequence nobody documented:
every acceleration the solver reasons about — including the `beta` setpoint
the user supplies — is the **gravity-offset pseudo-acceleration**
`Ẍ' = Ẍ_true - a_g`, not the true base-frame acceleration.

The class documentation says the opposite thing to a casual reader. It states
`beta = alpha^T * X_dotdot_N` and calls `X_dotdot_N` "the desired
(task-defined) spatial acceleration", with no mention of a gravity offset. A
user who wants the end-effector to *hold still* naturally writes `beta = 0`
and gets free-fall instead. This has already caused a downstream fork to patch
the solver's `constraint_calculation` to change the convention (see
"Maintenance notes"), when a one-line change at the call site would have done.

This plan does two things: it removes the fork's `E_input` patch so the solver
matches upstream Eq. (3.37) again (Step 1), and it documents the actual
convention with a test that pins it (Steps 2–3), so the next person does not
have to rediscover it with a debugger. The caller-side compensation that
replaces `E_input` lands in the `motion-spec` template PR, not here.

## Current state

Files:

- `orocos_kdl/src/chainhdsolver_vereshchagin.hpp` — the solver's class
  documentation lives entirely in this header (a deliberate upstream choice,
  see commit `b9570d6` "Documentation updated and moved to the header file").
  The `beta` description is at lines 192–214.
- `orocos_kdl/src/chainhdsolver_vereshchagin.cpp` — `constraint_calculation()`
  is where `beta` enters the balance. Step 1 removes the two fork-local
  `E_input` lines here; nothing else in this file changes.
  **Line numbers below are upstream's.** On the fork branch every line after
  the `E_input = s.E_tilde;` assignment is shifted by +1. Match on code
  content, not line number.
- `orocos_kdl/tests/solvertest.cpp` — `VereshchaginTest()` at line 818.

The documentation as it exists today, `chainhdsolver_vereshchagin.hpp:192-214`:

```
 * The last example involves the full specification of the desired end-effector motion (in this case,
 * not necessarily zero accelerations), i.e. specification of constraints in all 6 **DOFs**:
 *
 * **alpha** = ... (6x6 identity table) ...
 *
 * **beta** = **alpha^T * X_dotdot_N**
 *
 * Here, **N** stands for the index of the last robot's segment, end-effector (tool-tip). The reader should note
 * that we can directly assign values (magnitudes) of the desired (task-defined) spatial acceleration **6 x 1**
 * vector **X_dotdot_N** to the **6 x 1** vector of acceleration energy (**beta**) [3]. Even though physical
 * dimensions (units) of these two vectors are not the same, the property of matrix **alpha** (it contains
 * **unit** vectors), permits that we can assign values of desired accelerations to acceleration energy setpoints,
 * in respective directions.
```

The balance that consumes it, `chainhdsolver_vereshchagin.cpp:288-296`:

```cpp
    Vector6d acc;
    acc << Eigen::Vector3d::Map(acc_root.rot.data), Eigen::Vector3d::Map(acc_root.vel.data);
    nu_sum.noalias() = -(results[0].E_tilde.transpose() * acc);
    //nu_sum.setZero();
    nu_sum += beta.data;
    nu_sum -= results[0].G;

    //equation f) nu = M_0_inverse*(beta_N - E0_tilde`*acc0 - G0)
    nu.noalias() = M_0_inverse * nu_sum;
```

This matches Shakhimardanov 2015 Eq. (3.37) exactly: `L_0 ν = b_N - A_0^T Ẍ_0 - U_0`.
`Ẍ_0` here is `acc_root`, and every `Ẍ_i` downstream of it therefore carries
the `-a_g` offset. The constraint is imposed on the offset acceleration.

Measured evidence (3-link planar arm, `root_acc = (0,0,+9.81)`, constraints on
linear x, linear z, angular y; produced by `plans/reference/achd_beta_semantics.cpp`):

```
[beta = 0]                        qdd=[  41.495  -82.989   41.495]   <- end-effector free-falls
[beta = alpha^T * root_acc]       qdd=[  -0.000    0.000   -0.000]   <- end-effector holds still
```

The existing upstream test is self-consistent with this and must keep passing.
`solvertest.cpp:962-968` asserts `beta == getTransformedLinkAcceleration()`,
and `getTransformedLinkAcceleration()` also returns the offset acceleration, so
both sides carry the same offset:

```cpp
    std::vector<Twist> xDotdot(ns + 1);
    // This solver's function returns Cartesian accelerations of links in robot base coordinates
    vereshchaginSolver.getTransformedLinkAcceleration(xDotdot);
    CPPUNIT_ASSERT(Equal(beta_energy(0), xDotdot[ns].vel(0), eps));
    CPPUNIT_ASSERT(Equal(beta_energy(1), xDotdot[ns].vel(1), eps));
    CPPUNIT_ASSERT(Equal(beta_energy(2), xDotdot[ns].vel(2), eps));
    CPPUNIT_ASSERT(Equal(beta_energy(5), xDotdot[ns].rot(2), eps));
```

Repo conventions to match:

- Documentation is Doxygen block comments in the header, Markdown-flavoured,
  using `**bold**` for interface names and `[N]` bracketed reference numbers
  that resolve against the `## REFERENCES` list at the end of the same header.
  Match that style exactly — see the surrounding `#### External Forces: f_ext`
  section at `chainhdsolver_vereshchagin.hpp:216-222` as the exemplar.
- Tests are CppUnit. New assertions go inside an existing `void
  SolverTest::<Name>()` in `orocos_kdl/tests/solvertest.cpp`; new test methods
  must additionally be declared in `orocos_kdl/tests/solvertest.hpp` and
  registered with `CPPUNIT_TEST(<Name>);` there (see `solvertest.hpp:36`).

## Commands you will need

| Purpose | Command | Expected on success |
|---|---|---|
| Configure | `cd orocos_kdl && mkdir -p build && cd build && cmake -DENABLE_TESTS:BOOL=ON -DCMAKE_BUILD_TYPE=Release ./..` | exit 0 |
| Build | `cd orocos_kdl/build && make -j$(nproc)` | exit 0, no warnings from the files you touched |
| Test | `cd orocos_kdl/build && make check` | exit 0, all CppUnit tests pass |
| Single test | `cd orocos_kdl/build && ./tests/solvertest` | `OK (N tests)` |

Dependencies (already required by CI): `libeigen3-dev`, `libcppunit-dev`.

## Scope

**In scope** (the only files you may modify):
- `orocos_kdl/src/chainhdsolver_vereshchagin.cpp` — **Step 1 only**: delete the
  two `E_input` lines. No other change to this file.
- `orocos_kdl/src/chainhdsolver_vereshchagin.hpp` — the `E_input` member
  declaration (Step 1) and documentation comments (Step 2).
- `orocos_kdl/tests/solvertest.cpp`
- `orocos_kdl/tests/solvertest.hpp`

**Out of scope** (do NOT touch):
- The SVD truncation threshold `if (Sm(i) < 1e-8)`. Upstream has `1e-14`; the
  fork's value stays as-is by explicit decision. Leave it alone.
- Any other line of `constraint_calculation()` beyond the `E_input` removal.
  The balance must end up matching upstream exactly, not "improved".
- Any `*_fext*` or `*_fixed_joint*` solver file — these are fork-local files
  that do not exist upstream and are retired by plan 004.
- `getTransformedLinkAcceleration()` and the other getters — their offset
  convention is consistent with `beta` today and changing one without the
  other would desynchronise them.

## Git workflow

- Branch: `docs/vereshchagin-beta-gravity-convention`
- One commit is fine. Message style matches `git log --oneline` on this repo,
  e.g. `[Vereshchagin solver] Documentation updated and moved to the header file`.
  Suggested: `[Vereshchagin solver] Document the gravity-offset convention of beta`
- Do NOT push or open a PR unless the operator instructed it.

## Steps

### Step 1: Remove the `E_input` gravity shift (commit 1, fork-only)

Delete all three `E_input` sites so `constraint_calculation()` matches upstream
Eq. (3.37) again:

1. `orocos_kdl/src/chainhdsolver_vereshchagin.hpp` — the member declaration
   `Matrix6Xd E_input; //` in the private section.
2. `orocos_kdl/src/chainhdsolver_vereshchagin.cpp`, in `downwards_sweep()`'s
   `i == (int)ns` leaf branch — the line `E_input = s.E_tilde;`.
3. `orocos_kdl/src/chainhdsolver_vereshchagin.cpp`, in
   `constraint_calculation()` — the line `nu_sum += E_input.transpose() * acc;`.

Leave the SVD threshold alone.

After this, `constraint_calculation()`'s balance must read exactly:

```cpp
    nu_sum.noalias() = -(results[0].E_tilde.transpose() * acc);
    //nu_sum.setZero();
    nu_sum += beta.data;
    nu_sum -= results[0].G;
```

Commit this on its own, message e.g.
`[Vereshchagin solver] Remove non-upstream gravity shift from the constraint balance`.
This commit is **not** cherry-picked upstream — it only undoes a fork-local
change.

**Verify**: `grep -rn "E_input" orocos_kdl/src/` → no output.

**Verify**: `git diff origin/master -- orocos_kdl/src/chainhdsolver_vereshchagin.cpp`
→ the only remaining differences are the `1e-8` threshold and whatever plans
002/003 added. No `E_input`.

**Verify**: `cd orocos_kdl/build && make -j$(nproc)` → exit 0.

Note that `make check` may now FAIL on downstream expectations if anything was
tuned against the shifted convention — that is expected and is exactly why this
plan is round 2. `VereshchaginTest` itself must still pass, because it was
written against upstream's convention.

### Step 2: Add a gravity-offset subsection to the `beta` documentation

In `orocos_kdl/src/chainhdsolver_vereshchagin.hpp`, immediately after the
sentence ending `...in respective directions.` (currently line 214) and before
the `#### External Forces: f_ext` heading (currently line 216), insert a new
paragraph. It must state, in the header's existing Doxygen/Markdown style:

- The solver takes gravity into account by setting the root segment's
  acceleration to the **negative** of the gravitational acceleration (this is
  why the constructor doc already warns the sign is opposite to the FD and RNE
  solvers — cross-reference that existing note).
- As a consequence every Cartesian acceleration inside the solver, including
  the one constrained by **alpha**/**beta** and the one returned by
  **getTransformedLinkAcceleration**, is the *gravity-offset* acceleration
  `X_dotdot' = X_dotdot - a_g`, not the acceleration measured in an inertial
  base frame.
- Therefore, to constrain the end-effector's **true** base-frame acceleration
  to `X_dotdot_desired`, the caller supplies
  `beta = alpha^T * (X_dotdot_desired - a_g)`, which with the KDL sign
  convention is `beta = alpha^T * X_dotdot_desired + alpha^T * root_acc`.
- Give the concrete common case: to hold the end-effector still under gravity,
  `beta = alpha^T * root_acc`, **not** `beta = 0`. `beta = 0` commands
  free-fall.
- Reference `[4]` (Featherstone) for the root-acceleration technique, since
  that entry already exists in the header's reference list.

Do not reflow or reindent any surrounding lines — keep the diff to inserted
lines only.

**Verify**: `git diff --stat orocos_kdl/src/chainhdsolver_vereshchagin.hpp`
→ shows only insertions (`+N`), zero deletions (`-0`).

**Verify**: `cd orocos_kdl/build && make -j$(nproc)` → exit 0.

### Step 3: Add a test that pins the convention

In `orocos_kdl/tests/solvertest.cpp`, add a new test method
`void SolverTest::VereshchaginGravityOffsetTest()`. Model its structure on the
existing `SolverTest::VereshchaginTest()` at line 818 — same includes, same
`CPPUNIT_ASSERT(Equal(...))` idiom, same `eps` tolerance already in scope.

The test builds a small chain and asserts both halves of the convention:

1. Build a 3-segment chain of `Joint::RotY` segments, each with frame
   `Frame(Vector(0,0,0.4))` and a `RigidBodyInertia` of mass 2.0, COG
   `Vector(0,0,0.2)`, `RotationalInertia(0.02666667, 0.02666667, 1e-4)`.
2. `Twist root_acc(Vector(0.0, 0.0, 9.81), Vector::Zero());`
3. `nc = 3`; `alpha` columns: `Twist(Vector(1,0,0), Vector::Zero())`,
   `Twist(Vector(0,0,1), Vector::Zero())`, `Twist(Vector::Zero(), Vector(0,1,0))`.
4. Joint positions `q = {0.3, -0.6, 0.4}`, zero velocity, zero feed-forward
   torque, zero external wrench.
5. **Assertion A — `beta = 0` does NOT hold the arm still.** Solve with
   `beta = 0` and assert `q_dotdot(0)` is far from zero (use
   `CPPUNIT_ASSERT(std::abs(q_dotdot(0)) > 1.0)`). This is the trap the
   documentation now warns about; if it ever silently becomes false, the
   convention changed and the docs are stale.
6. **Assertion B — `beta = alpha^T * root_acc` DOES hold the arm still.**
   Solve with `beta = (0.0, 9.81, 0.0)` and assert each
   `q_dotdot(i)` equals `0.0` within `eps`.

Then declare `void VereshchaginGravityOffsetTest();` in
`orocos_kdl/tests/solvertest.hpp` alongside the existing declarations
(near line 57) and register it with `CPPUNIT_TEST(VereshchaginGravityOffsetTest );`
next to the existing `CPPUNIT_TEST(VereshchaginTest );` at line 36 — match the
file's existing spacing style, which puts a space before the closing paren.

**Verify**: `cd orocos_kdl/build && make -j$(nproc) && ./tests/solvertest`
→ `OK (N tests)` where N is one greater than before your change. Record the
before and after numbers.

### Step 4: Confirm no existing behaviour moved

**Verify**: `cd orocos_kdl/build && make check` → exit 0. In particular
`VereshchaginTest` and `FdAndVereshchaginSolversConsistencyTest` must still
pass with their existing hard-coded expected values (e.g. `nu(0) ==
669693.30355`). If either fails, you changed behaviour — STOP.

## Test plan

- New test: `VereshchaginGravityOffsetTest` in `orocos_kdl/tests/solvertest.cpp`,
  covering the two cases in Step 3 (`beta = 0` → not static; `beta = alpha^T *
  root_acc` → static).
- Structural pattern to follow: `SolverTest::VereshchaginTest()` at
  `orocos_kdl/tests/solvertest.cpp:818`.
- Verification: `cd orocos_kdl/build && make check` → all pass, including the
  one new test.

## Done criteria

ALL must hold:

- [ ] `cd orocos_kdl/build && make -j$(nproc)` exits 0
- [ ] `cd orocos_kdl/build && make check` exits 0
- [ ] `./tests/solvertest` reports exactly one more test than before the change
- [ ] `git diff orocos_kdl/src/chainhdsolver_vereshchagin.hpp` contains only
      added comment lines — verify with
      `git diff -U0 orocos_kdl/src/chainhdsolver_vereshchagin.hpp | grep '^-[^-]'`
      returning no output
- [ ] `git status --porcelain` lists only the three in-scope files
- [ ] `plans/README.md` status row updated

## STOP conditions

Stop and report back (do not improvise) if:

- The excerpts in "Current state" do not match the live code.
- `VereshchaginTest` or `FdAndVereshchaginSolversConsistencyTest` fails at any
  point. You are not authorised to change their expected values in this plan.
- Assertion A in Step 3 fails — i.e. `beta = 0` *does* hold the arm still.
  That would mean the convention is not what this plan documents. Report the
  measured `q_dotdot` values instead of adjusting the assertion.
- You conclude the `.cpp` must change to make the docs true. It must not; the
  docs are being made true to the code, not the reverse.

## Maintenance notes

- **This is the upstream half of a two-sided issue.** The downstream fork
  (`secorolab/orocos_kinematics_dynamics`, commit `42a9af9`) patched
  `constraint_calculation` to add `nu_sum += E_input.transpose() * acc`, which
  silently redefines `beta` to mean the true acceleration. That change was
  measured to be *exactly* equivalent to the caller passing
  `beta + alpha^T * root_acc` (max deviation `1.4e-14` in `tau` and `q_dotdot`
  across a 12-pose sweep with non-zero velocity, feed-forward torque and
  external wrench — see `plans/reference/achd_beta_semantics.cpp`). It is a
  change of input convention, not a correctness fix, and it breaks the
  `beta == getTransformedLinkAcceleration()` correspondence that
  `VereshchaginTest` asserts. Plan 004 removes it from the fork in favour of
  the caller-side shift this plan documents.
- A reviewer should check that the new paragraph does not contradict the
  constructor's existing note about the root-acceleration sign, and that the
  `[4]` reference number still points at Featherstone in the reference list.
- Deferred deliberately: making `getTransformedLinkAcceleration()` return the
  true acceleration, or adding an opt-in convention flag. Both are behaviour
  changes and need their own discussion with upstream maintainers.
