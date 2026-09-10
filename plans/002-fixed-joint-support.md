# Plan 002: Complete fixed-joint support in ChainHdSolver_Vereshchagin

> **Executor instructions**: Follow this plan step by step. Run every
> verification command and confirm the expected result before moving to the
> next step. If anything in the "STOP conditions" section occurs, stop and
> report — do not improvise. When done, update the status row for this plan
> in `plans/README.md`.
>
> **Drift check (run first)**:
> `git diff --stat 4cadf46..HEAD -- orocos_kdl/src/chainhdsolver_vereshchagin.cpp orocos_kdl/tests/solvertest.cpp`
> If any in-scope file changed since this plan was written, compare the
> "Current state" excerpts against the live code before proceeding; on a
> mismatch, treat it as a STOP condition.

## Status

- **Priority**: P1
- **Effort**: M
- **Risk**: MED
- **Depends on**: none (independent of 001, but land 001 first — it is trivial)
- **Category**: bug
- **Planned at**: commit `4cadf46`, 2026-08-07

## Why this matters

`ChainHdSolver_Vereshchagin` *appears* to support chains containing
`Joint::Fixed` segments — there are four `!= Joint::Fixed` guards in the file.
But all four only guard the joint-index counter `j`. The numerical paths are
unguarded, and a fixed joint makes them degenerate:

- `Z = 0` for a fixed joint, so `D = dot(Z, PZ) = 0`.
- Three separate places divide by `D`. Each produces `NaN`.
- The `NaN` enters `s.acc` and propagates to every segment after it, so the
  whole solve returns `NaN` — silently, with `E_NOERROR`.

There is also a latent out-of-bounds read: `s.F = segment.pose(q(j))` at
`chainhdsolver_vereshchagin.cpp:93` runs *before* the `j++` guard at line 132,
so on a fixed segment it reads `q(j)` for the wrong joint — and when the fixed
segment comes after the last movable joint, `j == nj` and it reads past the end
of the `JntArray`.

This matters in practice because URDF-derived chains routinely carry fixed
segments (tool flanges, sensor mounts, gripper adapters). Today the only way to
use this solver on such a chain is to pre-collapse the fixed joints yourself.

## Current state

Files:

- `orocos_kdl/src/chainhdsolver_vereshchagin.cpp` — all four degenerate sites.
- `orocos_kdl/tests/solvertest.cpp` — `VereshchaginTest()` at line 818 is the
  structural exemplar for the new test.
- `orocos_kdl/tests/solvertest.hpp` — test registration, line 36.

**Site 1 — `initial_upwards_sweep`, `chainhdsolver_vereshchagin.cpp:89-133`.**
`q(j)` is read three times before the guard at the bottom:

```cpp
        const Segment& segment = chain.getSegment(i);
        segment_info& s = results[i + 1];
        //The pose between the joint root and the segment tip (tip expressed in joint root coordinates)
        s.F = segment.pose(q(j)); //X pose of each link in link coord system

        F_total = F_total * s.F; //X pose of the each link in root coord system
        s.F_base = F_total; //X pose of the each link in root coord system for getter functions

        //The velocity due to the joint motion of the segment expressed in the segments reference frame (tip)
        Twist vj = s.F.M.Inverse(segment.twist(q(j), qdot(j))); //XDot of each link

        //The unit velocity due to the joint motion of the segment expressed in the segments reference frame (tip)
        s.Z = s.F.M.Inverse(segment.twist(q(j), 1.0));
        //Put Z in the joint root reference frame:
        s.Z = s.F * s.Z;
        ...
        Wrench FextLocal = F_total.M.Inverse() * f_ext[i];
        s.U = s.v * (s.H * s.v) - FextLocal; //f_ext[i];
        if (segment.getJoint().getType() != Joint::Fixed)
            j++;
```

**Site 2 — `downwards_sweep`, `chainhdsolver_vereshchagin.cpp:178-213`.** Four
divisions by `child.D`:

```cpp
            segment_info& child = results[i + 1];
            //Copy PZ into a vector so we can do matrix manipulations, put torques above forces
            Vector6d vPZ;
            vPZ << Eigen::Vector3d::Map(child.PZ.torque.data), Eigen::Vector3d::Map(child.PZ.force.data);
            Matrix6d PZDPZt;
            PZDPZt.noalias() = vPZ * vPZ.transpose();
            PZDPZt /= child.D;

            s.P_tilde = s.H + child.P - ArticulatedBodyInertia(PZDPZt.bottomRightCorner<3,3>(), PZDPZt.topRightCorner<3,3>(), PZDPZt.topLeftCorner<3,3>());
            s.R_tilde = s.U + child.R + child.PC + (child.PZ / child.D) * child.u;
            s.E_tilde = child.E;
            s.E_tilde.noalias() -= (vPZ * child.EZ.transpose()) / child.D;

            s.M = child.M;
            s.M.noalias() -= (child.EZ * child.EZ.transpose()) / child.D;

            s.G = child.G;
            Twist CiZDu = child.C + (child.Z / child.D) * child.u;
            Vector6d vCiZDu;
            vCiZDu << Eigen::Vector3d::Map(CiZDu.rot.data), Eigen::Vector3d::Map(CiZDu.vel.data);
            s.G.noalias() += child.E.transpose() * vCiZDu;
```

**Site 3 — `downwards_sweep`, `chainhdsolver_vereshchagin.cpp:238-259`.**
`ff_torques(j)` is read without a fixed-joint guard:

```cpp
            if (chain.getSegment(i - 1).getJoint().getType() != Joint::Fixed)
                s.D = chain.getSegment(i - 1).getJoint().getInertia() + dot(s.Z, s.PZ);
            else
                s.D = dot(s.Z, s.PZ);

            s.PC = s.P * s.C;

            //projection of coriolis and centrepital forces into joint subspace (0 0 Z)
            s.totalBias = -dot(s.Z, s.R + s.PC);
            s.u = ff_torques(j) + s.totalBias;
            ...
            if (chain.getSegment(i - 1).getJoint().getType() != Joint::Fixed)
                j--;
```

**Site 4 — `final_upwards_sweep`, `chainhdsolver_vereshchagin.cpp:320-347`.**
Three divisions by `s.D` and three writes indexed by `j`, all before the guard:

```cpp
        Wrench parent_force = s.P*a_p;
        double parent_forceProjection = -dot(s.Z, parent_force);
        double parentAccComp = parent_forceProjection / s.D;

        constraint_torques(j) = -dot(s.Z, constraint_force);
        total_torques(j) = s.u + parent_forceProjection + constraint_torques(j);

        s.constAccComp = constraint_torques(j) / s.D;
        s.nullspaceAccComp = s.u / s.D;

        q_dotdot(j) = (s.nullspaceAccComp + parentAccComp + s.constAccComp);
        s.acc = s.F.Inverse(a_p + s.Z * q_dotdot(j) + s.C);
        if (chain.getSegment(i - 1).getJoint().getType() != Joint::Fixed)
            j++;
```

**The correct limit.** A fixed joint contributes no joint-space degree of
freedom. Every term above that divides by `D` is a projection onto the joint
axis, and vanishes when `Z = 0`. Physically this is a rigid-body merge: the
child's articulated inertia and bias force add to the parent's with no
reduction. Concretely:

| Quantity | Movable joint | Fixed joint |
|---|---|---|
| `P_tilde` (Eq. 3.31) | `H + child.P - PZ·PZ^T/D` | `H + child.P` |
| `R_tilde` (Eq. 3.32) | `U + child.R + child.PC + (child.PZ/D)·child.u` | `U + child.R + child.PC` |
| `E_tilde` (Eq. 3.33) | `child.E - vPZ·child.EZ^T/D` | `child.E` |
| `M` (Eq. 3.35) | `child.M - child.EZ·child.EZ^T/D` | `child.M` |
| `G` (Eq. 3.34) | `child.G + child.E^T·(child.C + child.Z/D·child.u)` | `child.G + child.E^T·child.C` |
| `q_dotdot`, `acc` (Eq. 3.30) | `D^-1{...}`; `acc = F^-1(a_p + Z·q̈ + C)` | no `q̈`; `acc = F^-1(a_p + C)` |

Equation numbers refer to Shakhimardanov, "Composable robot motion stack",
PhD thesis, KU Leuven 2015 — already reference `[3]` in the solver header.

Repo conventions: 4-space indent, `//` comments in the existing terse style,
CppUnit tests (see plan 001's "Current state" for the registration mechanics).

## Commands you will need

| Purpose | Command | Expected on success |
|---|---|---|
| Configure | `cd orocos_kdl && mkdir -p build && cd build && cmake -DENABLE_TESTS:BOOL=ON -DCMAKE_BUILD_TYPE=Release ./..` | exit 0 |
| Build | `cd orocos_kdl/build && make -j$(nproc)` | exit 0 |
| Test | `cd orocos_kdl/build && make check` | exit 0, all pass |
| Single test | `cd orocos_kdl/build && ./tests/solvertest` | `OK (N tests)` |

## Scope

**In scope**:
- `orocos_kdl/src/chainhdsolver_vereshchagin.cpp`
- `orocos_kdl/tests/solvertest.cpp`
- `orocos_kdl/tests/solvertest.hpp`

**Out of scope** (do NOT touch):
- `orocos_kdl/src/chainhdsolver_vereshchagin.hpp` — no API change is needed
  or wanted here. The fix is entirely internal.
- The `nc == 0` / unconstrained code path — orthogonal, do not "improve" it.
- `orocos_kdl/src/joint.hpp`, `orocos_kdl/src/segment.hpp` — the downstream
  fork added `setInertia()` and `getMutableJoint()` accessors to these. They
  are **not** needed for this plan and must not be pulled in.
- Any `*_fixed_joint*` file — those are fork-local duplicates. This plan fixes
  the real solver in place instead of adding a fifth copy.

## Git workflow

- Branch: `fix/vereshchagin-fixed-joint-support`
- Commit per site or one commit for the fix plus one for tests; message style
  matches `git log --oneline`, e.g.
  `[Vereshchagin solver] Handle fixed joints without producing NaN`
- Do NOT push or open a PR unless the operator instructed it.

## Steps

### Step 1: Add a failing test first

Before changing any solver code, add `void SolverTest::VereshchaginFixedJointTest()`
to `orocos_kdl/tests/solvertest.cpp`, declared in `solvertest.hpp` and
registered with `CPPUNIT_TEST(VereshchaginFixedJointTest );` next to line 36.
Follow the structure of `VereshchaginTest()` at line 818.

The test builds **two chains that are physically identical**:

- Chain A: 3 segments, all `Joint::RotY`, frame `Frame(Vector(0,0,0.4))`,
  `RigidBodyInertia` mass 2.0, COG `Vector(0,0,0.2)`,
  `RotationalInertia(0.02666667, 0.02666667, 1e-4)`.
- Chain B: the same, but with a **massless fixed segment inserted between the
  2nd and 3rd movable segments** — `Segment(Joint(Joint::Fixed),
  Frame(Vector(0,0,0.0)), RigidBodyInertia::Zero())`. A zero-length, zero-mass
  fixed segment changes nothing physically, so both chains must produce
  identical results.

Solve both with the same `root_acc = Twist(Vector(0,0,9.81), Vector::Zero())`,
`nc = 3`, alpha columns `Twist(Vector(1,0,0),Vector::Zero())`,
`Twist(Vector(0,0,1),Vector::Zero())`, `Twist(Vector::Zero(),Vector(0,1,0))`,
non-zero `q`, `qdot`, `ff_torques` and a non-zero `f_ext` on the last segment.

Assert, for all 3 joints:
- `q_dotdot_A(i) == q_dotdot_B(i)` within `eps`
- `constraint_torques_A(i) == constraint_torques_B(i)` within `eps`
- and explicitly that no result is NaN:
  `CPPUNIT_ASSERT(q_dotdot_B(i) == q_dotdot_B(i));`

Note the array sizes: chain B has 4 segments but still 3 joints, so
`q`, `qdot`, `q_dotdot`, `ff_torques`, `constraint_torques` are all size 3
while `f_ext` is size 4.

**Verify**: `cd orocos_kdl/build && make -j$(nproc) && ./tests/solvertest`
→ the new test **FAILS** (expect NaN mismatches on chain B). Record the
failure output. If it passes, STOP — the premise of this plan is wrong.

### Step 2: Guard `initial_upwards_sweep`

In `chainhdsolver_vereshchagin.cpp`, hoist a
`const bool segment_is_fixed = (segment.getJoint().getType() == Joint::Fixed);`
immediately after `segment_info& s = results[i + 1];`, then branch:

- Movable: keep the existing body verbatim, and move the `j++` into the end of
  this branch (removing the trailing `if (...) j++;`).
- Fixed: `s.F = segment.pose(0.0);` (KDL ignores the argument for a fixed
  joint — this exists only to avoid the out-of-bounds `q(j)` read), then
  `F_total`/`s.F_base` as before, `Twist vj = Twist::Zero();`,
  `s.Z = Twist::Zero();`, and the velocity recursion without the joint term:
  `s.v = s.F.Inverse(results[i].v)` for `i != 0`, else `Twist::Zero()`.
  `s.A` is unchanged in both branches. Do **not** increment `j`.

`vj` must be declared before the branch so the `s.C = s.v * vj;` line below
still compiles and yields `Twist::Zero()` for a fixed joint.

Leave `s.C`, `s.H` and `s.U` computation shared below the branch — they are
correct for both cases.

**Verify**: `cd orocos_kdl/build && make -j$(nproc)` → exit 0.

### Step 3: Guard the child recursion in `downwards_sweep`

Inside the `else` branch at line 178, before the `vPZ` computation, add:

```cpp
const bool child_is_fixed = (chain.getSegment(i).getJoint().getType() == Joint::Fixed);
```

Note the index: `child` is `results[i + 1]`, whose joint is
`chain.getSegment(i)`. Getting this index wrong is the most likely bug in this
plan — double-check it against the existing `results[i + 1]` line.

Branch: keep the existing five assignments for the movable case; for the fixed
case use the right-hand column of the table in "Current state" — no division by
`child.D` anywhere.

**Verify**: `grep -n "child.D" orocos_kdl/src/chainhdsolver_vereshchagin.cpp`
→ every hit is inside the `if (!child_is_fixed)` branch.

### Step 4: Guard `ff_torques(j)` in `downwards_sweep`

At line 250-251, replace the unconditional `s.u = ff_torques(j) + s.totalBias;`
with a branch on `chain.getSegment(i - 1).getJoint().getType() != Joint::Fixed`.
For a fixed joint, `s.u = s.totalBias;` — which is `0`, since `Z = 0` makes
`totalBias = -dot(Z, R + PC) = 0`. The point of the branch is to avoid indexing
`ff_torques(j)` with a `j` belonging to a different joint.

**Verify**: `cd orocos_kdl/build && make -j$(nproc)` → exit 0.

### Step 5: Guard `final_upwards_sweep`

Wrap the body from `Vector6d tmp = s.E*nu;` through `s.acc = ...` in
`if (chain.getSegment(i - 1).getJoint().getType() != Joint::Fixed) { ... j++; }`,
moving the existing `j++` inside. In the `else` branch:

```cpp
s.constAccComp = 0.0;
s.nullspaceAccComp = 0.0;
s.acc = s.F.Inverse(a_p + s.C);
```

No writes to `q_dotdot`, `constraint_torques` or `total_torques`, and no `j++`.

Also delete the unused `Twist a_g;` declaration at line 309 while you are in
this function — it is dead and the compiler warns on it.

**Verify**: `cd orocos_kdl/build && make -j$(nproc) && ./tests/solvertest`
→ `OK (N tests)`, including the new `VereshchaginFixedJointTest`.

### Step 6: Confirm no regression

**Verify**: `cd orocos_kdl/build && make check` → exit 0. `VereshchaginTest`
must still produce its hard-coded expected values (`nu(0) == 669693.30355`
etc.) — those chains have no fixed joints, so every new branch must be
untaken and the numbers must be bit-identical.

## Test plan

- New test `VereshchaginFixedJointTest` in `orocos_kdl/tests/solvertest.cpp`:
  equivalence of a chain with an inserted null fixed segment against the same
  chain without it, on `q_dotdot` and `constraint_torques`, plus explicit
  NaN checks.
- Structural pattern: `SolverTest::VereshchaginTest()` at
  `orocos_kdl/tests/solvertest.cpp:818`.
- Verification: `cd orocos_kdl/build && make check` → all pass.

## Done criteria

ALL must hold:

- [ ] `cd orocos_kdl/build && make check` exits 0
- [ ] `VereshchaginFixedJointTest` exists and passes
- [ ] `VereshchaginTest` and `FdAndVereshchaginSolversConsistencyTest` pass
      with unchanged expected values
- [ ] `grep -c "Joint::Fixed" orocos_kdl/src/chainhdsolver_vereshchagin.cpp`
      returns at least 6 (was 4)
- [ ] No new file was created under `orocos_kdl/src/`
      (`git status --porcelain orocos_kdl/src/` shows only the modified `.cpp`)
- [ ] `plans/README.md` status row updated

## STOP conditions

Stop and report back (do not improvise) if:

- The Step 1 test **passes** before any solver change. The premise (fixed
  joints currently produce NaN) would be false.
- After Step 5 the two chains still disagree by more than `eps`. Do not widen
  the tolerance to make it pass — report the actual deviation per joint.
- `VereshchaginTest`'s hard-coded values change at all. That means a branch is
  being taken on a chain with no fixed joints, which is a logic error in your
  guard conditions.
- You find you need `getMutableJoint()` or `Joint::setInertia()`. You do not;
  those are out of scope and their absence is not a blocker.

## Maintenance notes

- The downstream fork already carries this fix, as two duplicated solver files
  (`chainhdsolver_vereshchagin_fixed_joint.{hpp,cpp}` and the `_fext_` pair),
  totalling ~2100 lines of copy-paste. Landing this plan lets plan 004 delete
  all of them. When comparing, note the fork's copies also contain unrelated
  changes (an `E_input` term and an SVD tolerance change) that must **not**
  come along — see plan 001.
- A reviewer should scrutinise the index arithmetic in Step 3
  (`chain.getSegment(i)` for `results[i + 1]`) and confirm that `j` is
  incremented in exactly one place per sweep.
- Deferred: `Joint::None` is a distinct enum value from `Joint::Fixed` in KDL
  and is *not* handled by this plan. If a chain uses it, behaviour is
  unchanged from today. Worth a follow-up issue but not this PR.
