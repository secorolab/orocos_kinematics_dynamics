# Plan 003: Add per-driver acceleration-energy weighting to ChainHdSolver_Vereshchagin

> **Executor instructions**: Follow this plan step by step. Run every
> verification command and confirm the expected result before moving to the
> next step. If anything in the "STOP conditions" section occurs, stop and
> report — do not improvise. When done, update the status row for this plan
> in `plans/README.md`.
>
> **Drift check (run first)**:
> `git diff --stat 4cadf46..HEAD -- orocos_kdl/src/chainhdsolver_vereshchagin.cpp orocos_kdl/src/chainhdsolver_vereshchagin.hpp orocos_kdl/tests/solvertest.cpp`
> If any in-scope file changed since this plan was written, compare the
> "Current state" excerpts against the live code before proceeding; on a
> mismatch, treat it as a STOP condition.

## Status

- **Priority**: P1
- **Effort**: M
- **Risk**: MED
- **Depends on**: plans/002-fixed-joint-support.md (only to avoid conflicting
  edits in the same functions — the feature itself is independent)
- **Category**: direction
- **Planned at**: commit `4cadf46`, 2026-08-07

## Why this matters

`ChainHdSolver_Vereshchagin` hard-codes one prioritisation policy: Cartesian
acceleration constraints always win over external wrenches and feed-forward
joint torques. The solver's own header acknowledges this and says the
prioritisation "can be changed (see [3] & [5] for more details) but those
features are not implemented in KDL".

The consequence for users is that the solver actively *fights* an external
wrench. Measured on a 3-link arm: adding a 10 N wrench at the tip changes the
constraint torque by `(-1.60, +1.42, -0.40)` Nm — the solver compensates it
exactly. Any application that wants compliance (admittance control, contact
tasks, a deliberate support force) cannot express it, and users work around it
by running two solver instances and summing the torques. That workaround is
wrong: the second pass is contaminated by gravity and by velocity-product bias
terms, so the sum double-counts them.

Shakhimardanov 2015 §3.5.3 (reference `[3]` in the solver header) gives the
fix, and it is small. Weighting is applied to the acceleration-energy balance
that already exists in `constraint_calculation()`, Eq. (3.42):

> `b_control = w_posture * b_posture + w_ee * b_ee`

This plan implements exactly that, additively, with defaults that reproduce
today's arithmetic **line for line** — not merely to within a tolerance.

## Current state

Files:

- `orocos_kdl/src/chainhdsolver_vereshchagin.hpp` — class declaration; the
  `segment_info` struct is at lines 486–529, private members at 466–484.
- `orocos_kdl/src/chainhdsolver_vereshchagin.cpp` — the three sweeps.
- `orocos_kdl/tests/solvertest.cpp` / `.hpp` — CppUnit tests.

**The key structural fact this plan exploits.** In the inward recursion, the
*inertia* quantities are completely independent of any force input. Only
`R`/`u`/`G` depend on the drivers, and they are **affine** in them with
driver-independent coefficients:

| Quantity | Depends on drivers? |
|---|---|
| `P_tilde`, `P`, `PZ`, `D`, `Z`, `C`, `PC`, `E`, `E_tilde`, `EZ`, `M` | No |
| `U`, `R_tilde`, `R`, `u`, `G` | Yes, affinely |

So a per-driver contribution can be tracked through the *same* recursion in a
parallel accumulator, and the accumulators sum back to today's totals exactly.

**Where the drivers enter.** `chainhdsolver_vereshchagin.cpp:128-131`
(external wrench into `s.U`):

```cpp
        //wrench of the rigid body bias forces and the external forces on the segment (in body coordinates, tip)
        //external forces are taken into account through s.U.
        Wrench FextLocal = F_total.M.Inverse() * f_ext[i];
        s.U = s.v * (s.H * s.v) - FextLocal; //f_ext[i];
```

`chainhdsolver_vereshchagin.cpp:195` and `:207-212` (propagation into `R` and
the acceleration-energy accumulator `G`):

```cpp
            s.R_tilde = s.U + child.R + child.PC + (child.PZ / child.D) * child.u;
            ...
            //equation e) (see Vereshchagin89)
            s.G = child.G;
            Twist CiZDu = child.C + (child.Z / child.D) * child.u;
            Vector6d vCiZDu;
            vCiZDu << Eigen::Vector3d::Map(CiZDu.rot.data), Eigen::Vector3d::Map(CiZDu.vel.data);
            s.G.noalias() += child.E.transpose() * vCiZDu;
```

`chainhdsolver_vereshchagin.cpp:249-251` (feed-forward torque into `u`):

```cpp
            //projection of coriolis and centrepital forces into joint subspace (0 0 Z)
            s.totalBias = -dot(s.Z, s.R + s.PC);
            s.u = ff_torques(j) + s.totalBias;
```

**Where the weight applies** — `chainhdsolver_vereshchagin.cpp:288-296`. `G` is
the accumulated acceleration energy the non-constraint drivers have *already*
generated at the end-effector; it is subtracted from the demand `beta`, which
is precisely why the constraint compensates those drivers:

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

`segment_info` today, `chainhdsolver_vereshchagin.hpp:486-514` (abbreviated to
the members this plan touches):

```cpp
    struct segment_info
    {
        ...
        Wrench U; //wrench p of the bias forces (in cartesian space)
        Wrench R; //wrench p of the bias forces
        Wrench R_tilde; //vector of wrench p of the bias forces (new) in matrix form
        ...
        Eigen::VectorXd G; //magnitude of the constraint forces already generated at link i
        ...
        double u; //vector u[i] = torques(i) - S[i]^T*(p_A[i] + I_A[i]*C[i]) in joint subspace.

        segment_info(unsigned int nc):
            D(0),nullspaceAccComp(0),constAccComp(0),biasAccComp(0),totalBias(0),u(0)
        {
            E.resize(6, nc);
            E_tilde.resize(6, nc);
            G.resize(nc);
            M.resize(nc, nc);
            EZ.resize(nc);
            E.setZero();
            E_tilde.setZero();
            M.setZero();
            G.setZero();
            EZ.setZero();
        };
    };
```

Error-code convention, `chainhdsolver_vereshchagin.cpp:56-65`: size problems
return `(error = E_SIZE_MISMATCH)`. The enum lives in
`orocos_kdl/src/solveri.hpp` (`E_NOERROR = 0`, `E_SIZE_MISMATCH = -4`).

Repo conventions: 4-space indent; `Eigen::` fully qualified (upstream removed
`using namespace Eigen` in commit `86fffc1` — do not reintroduce it); public
API documented with Doxygen `/** ... */` blocks in the header.

## Design to implement

Track two extra accumulators alongside the existing totals — **the existing
`U`, `R`, `R_tilde`, `u`, `G` keep their present meaning and values unchanged**
(the full, unweighted totals), so `final_upwards_sweep()` and every getter need
no edit at all.

Two new channels per segment:

- `fext` channel: seeded with `U_fext = -FextLocal` (the external wrench only).
- `ff` channel: seeded with `U_ff = 0`; its `u_ff` picks up `ff_torques(j)`.

Each channel runs the same recursion, minus the terms that belong to the
"nature" part (the rigid-body bias `v x (H v)` and the velocity-product
acceleration `C`, which are physics, not drivers, and are never weighted):

```
R_tilde_c = U_c + child.R_c + (child.PZ / child.D) * child.u_c
R_c       = F * R_tilde_c
u_fext    =                 -dot(Z, R_fext)
u_ff      = ff_torques(j)   -dot(Z, R_ff)
G_c       = child.G_c + child.E^T * vec( (child.Z / child.D) * child.u_c )
```

Then in `constraint_calculation()`, credit only the weighted fraction of each
driver's already-generated acceleration energy, by *removing* the un-credited
part from the existing total:

```
nu_sum -= results[0].G
        - (1 - w_fext) ⊙ results[0].G_fext
        - (1 - w_ff)   ⊙ results[0].G_ff
```

With `w_fext = w_ff = 1` (the defaults) both correction terms are identically
zero and the executed arithmetic is byte-for-byte what it is today. That is why
the regression test can demand *exact* equality, not `eps` equality.

Semantics of the weights, for the header docs:

- `w = 1` — the constraint fully compensates that driver (today's behaviour).
- `w = 0` — the constraint is blind to that driver; the driver's effect passes
  through to the end-effector. This is the "prioritise the driver over the
  constraint" case of §3.5.3.
- `0 < w < 1` — partial credit, Eq. (3.42).

The weights are `nc`-vectors, not scalars: `beta` and `G` are both `nc x 1`
over the columns of `alpha`, so per-constraint-direction weighting is free and
is the useful case (be compliant along the contact normal, stay stiff in the
other directions).

## Commands you will need

| Purpose | Command | Expected on success |
|---|---|---|
| Configure | `cd orocos_kdl && mkdir -p build && cd build && cmake -DENABLE_TESTS:BOOL=ON -DCMAKE_BUILD_TYPE=Release ./..` | exit 0 |
| Build | `cd orocos_kdl/build && make -j$(nproc)` | exit 0 |
| Test | `cd orocos_kdl/build && make check` | exit 0, all pass |
| Single test | `cd orocos_kdl/build && ./tests/solvertest` | `OK (N tests)` |

## Scope

**In scope**:
- `orocos_kdl/src/chainhdsolver_vereshchagin.hpp`
- `orocos_kdl/src/chainhdsolver_vereshchagin.cpp`
- `orocos_kdl/tests/solvertest.cpp`
- `orocos_kdl/tests/solvertest.hpp`

**Out of scope** (do NOT touch):
- `final_upwards_sweep()` — by design this plan requires no change there.
  If you find yourself editing it, your channel bookkeeping is wrong. STOP.
- The `CartToJnt` signature — adding a parameter would break every caller and
  the PyKDL binding. Weights are set through a separate setter.
- `orocos_kdl/src/chainhdsolver_vereshchagin_fext*.{hpp,cpp}` and
  `*_fixed_joint*` — fork-local duplicates, retired by plan 004.
- `python_orocos_kdl/` — binding the new setter is a deliberate follow-up, not
  part of this PR.
- The SVD truncation policy and its `1e-14` threshold — unrelated.

## Git workflow

- Branch: `feature/vereshchagin-driver-weighting`
- Suggested commits: (1) channel accumulators, no behaviour change;
  (2) weights + setter; (3) tests. Message style matches `git log --oneline`,
  e.g. `[Vereshchagin solver] Add per-driver acceleration energy weighting`
- Do NOT push or open a PR unless the operator instructed it.

## Steps

### Step 1: Add the channel accumulators to `segment_info`

In `orocos_kdl/src/chainhdsolver_vereshchagin.hpp`, extend `segment_info` with:

```cpp
        // Per-driver contributions, tracked alongside the totals above so that
        // constraint_calculation() can credit each driver only partially.
        // Sum of the channels plus the rigid-body bias terms == the totals.
        Wrench U_fext, R_fext, R_tilde_fext;
        Wrench U_ff, R_ff, R_tilde_ff;
        double u_fext, u_ff;
        Eigen::VectorXd G_fext, G_ff;
```

Initialise `u_fext(0), u_ff(0)` in the constructor initialiser list next to the
existing `u(0)`, and inside the body resize and zero `G_fext` and `G_ff`
exactly as `G` is handled.

**Verify**: `cd orocos_kdl/build && make -j$(nproc)` → exit 0.

### Step 2: Seed the channels in `initial_upwards_sweep`

At `chainhdsolver_vereshchagin.cpp:128-131`, leave `s.U` exactly as it is and
add the channel seeds:

```cpp
        Wrench FextLocal = F_total.M.Inverse() * f_ext[i];
        s.U = s.v * (s.H * s.v) - FextLocal; //f_ext[i];
        s.U_fext = -FextLocal;      // external-wrench channel
        s.U_ff = Wrench::Zero();    // feed-forward channel is seeded at the joint, not the segment
```

**Verify**: `cd orocos_kdl/build && make -j$(nproc) && make check` → exit 0,
all tests still pass (nothing consumes the new members yet).

### Step 3: Propagate the channels in `downwards_sweep`

Three places, all mirroring the existing total:

1. In the `i == (int)ns` branch (leaf, around line 168–177), alongside
   `s.R_tilde = s.U;` add `s.R_tilde_fext = s.U_fext;`,
   `s.R_tilde_ff = s.U_ff;`, and zero `s.G_fext` / `s.G_ff` alongside the
   existing `s.G.setZero();`.

2. In the `else` branch, alongside line 195 and lines 207–212:

```cpp
            s.R_tilde_fext = s.U_fext + child.R_fext + (child.PZ / child.D) * child.u_fext;
            s.R_tilde_ff   = s.U_ff   + child.R_ff   + (child.PZ / child.D) * child.u_ff;
```

and, for the acceleration energy — note there is **no** `child.C` term here,
that belongs to the unweighted nature part:

```cpp
            s.G_fext = child.G_fext;
            Twist ZDu_fext = (child.Z / child.D) * child.u_fext;
            Vector6d vZDu_fext;
            vZDu_fext << Eigen::Vector3d::Map(ZDu_fext.rot.data), Eigen::Vector3d::Map(ZDu_fext.vel.data);
            s.G_fext.noalias() += child.E.transpose() * vZDu_fext;
```

and the same shape for `_ff`.

3. In the `if (i != 0)` block, alongside `s.R = s.F * s.R_tilde;` (around line
   218) add `s.R_fext = s.F * s.R_tilde_fext;` and
   `s.R_ff = s.F * s.R_tilde_ff;`. Then alongside lines 249–251, add the
   per-channel joint projections, leaving `s.u` untouched:

```cpp
            s.u_fext = -dot(s.Z, s.R_fext);
            s.u_ff   = ff_torques(j) - dot(s.Z, s.R_ff);
```

If plan 002 has landed, `s.u_ff` must sit inside the same non-fixed-joint
branch that guards `ff_torques(j)`, with `s.u_ff = 0.0;` in the fixed branch.

**Verify**: `cd orocos_kdl/build && make -j$(nproc) && make check` → exit 0.
Still no behaviour change: nothing reads the channels yet, so `VereshchaginTest`
must produce its exact existing values.

### Step 4: Add the weights and the public setter

In the header's public section, after `updateInternalDataStructures()`:

```cpp
    /**
     * Set the per-driver weights used when solving for the constraint force
     * magnitudes. Each weight is an nc-vector, one entry per column of alpha.
     *
     * A weight of 1.0 (the default) means the acceleration constraint fully
     * compensates that driver, i.e. the constraint is satisfied exactly
     * regardless of the driver -- this is the classic Popov-Vereshchagin
     * prioritisation. A weight of 0.0 means the constraint ignores that
     * driver, letting its effect pass through to the constrained segment.
     * Values in between blend the two, as in Eq. (3.42) of [3].
     *
     * \param w_f_ext weight per constraint direction for the external wrenches
     * \param w_ff_torques weight per constraint direction for the feed-forward joint torques
     * \return E_NOERROR on success, E_SIZE_MISMATCH if either vector is not of size nc
     */
    int setDriverWeights(const Eigen::VectorXd& w_f_ext, const Eigen::VectorXd& w_ff_torques);
```

Add private members `Eigen::VectorXd w_f_ext, w_ff_torques;`, initialise both
to `Eigen::VectorXd::Ones(nc)` in the constructor, and re-initialise them in
`updateInternalDataStructures()` (which currently only resizes
`total_torques` and `results` — match its style).

Implement the setter to validate `.size() == nc` for both, returning
`(error = E_SIZE_MISMATCH)` on failure and `(error = E_NOERROR)` otherwise.

**Verify**: `cd orocos_kdl/build && make -j$(nproc) && make check` → exit 0,
unchanged results.

### Step 5: Apply the weights in `constraint_calculation`

Replace the single line `nu_sum -= results[0].G;` at
`chainhdsolver_vereshchagin.cpp:293` with the total minus the un-credited
fractions:

```cpp
    // Credit each driver's already-generated acceleration energy only in
    // proportion to its weight. With unit weights both correction terms are
    // exactly zero and this reduces to nu_sum -= results[0].G.
    nu_sum -= results[0].G;
    nu_sum += (Eigen::VectorXd::Ones(nc) - w_f_ext).cwiseProduct(results[0].G_fext);
    nu_sum += (Eigen::VectorXd::Ones(nc) - w_ff_torques).cwiseProduct(results[0].G_ff);
```

Update the two `//equation f)` comments in this function to mention the
weighting, and cite `[3]` §3.5.3 / Eq. (3.42).

**Verify**: `cd orocos_kdl/build && make -j$(nproc) && make check` → exit 0.
`VereshchaginTest`'s hard-coded values (`nu(0) == 669693.30355` etc.) must be
**unchanged**, because the default weights make the added terms exactly zero.

### Step 6: Document the feature in the class header

In `chainhdsolver_vereshchagin.hpp`, the section
`### Prioritizations between motion drivers (interfaces)` (lines 277–292)
currently ends with:

```
 * Nevertheless, the above-described prioritization can be changed (see [3] & [5] for more details) but those features are not implemented in KDL.
```

Replace that sentence with a short subsection describing `setDriverWeights`:
the three weight regimes from the "Design to implement" section above, the fact
that weights are per constraint direction, that the default is all-ones and
reproduces the classic behaviour, and a pointer to `[3]` §3.5.3.

**Verify**: `grep -n "not implemented in KDL" orocos_kdl/src/chainhdsolver_vereshchagin.hpp`
→ no output.

### Step 7: Tests

Add `void SolverTest::VereshchaginDriverWeightingTest()` to
`orocos_kdl/tests/solvertest.cpp`, declared in `solvertest.hpp` and registered
with `CPPUNIT_TEST(VereshchaginDriverWeightingTest );` next to line 36. Model
on `SolverTest::VereshchaginTest()` at line 818. Use a 3-link `Joint::RotY`
chain (segment frame `Frame(Vector(0,0,0.4))`, `RigidBodyInertia` mass 2.0,
COG `Vector(0,0,0.2)`, `RotationalInertia(0.02666667, 0.02666667, 1e-4)`),
`root_acc = Twist(Vector(0,0,9.81), Vector::Zero())`, `nc = 3` with alpha
columns `Twist(Vector(1,0,0),Vector::Zero())`,
`Twist(Vector(0,0,1),Vector::Zero())`, `Twist(Vector::Zero(),Vector(0,1,0))`,
non-zero `q`, `qdot`, `ff_torques`, and a non-zero `f_ext` on the last segment.

Cases:

1. **Default weights are a no-op.** Solve once without calling
   `setDriverWeights`, once after calling it with two all-ones vectors. Assert
   the two `constraint_torques` and `q_dotdot` are **exactly** equal
   (`CPPUNIT_ASSERT_EQUAL(a, b)` on the doubles, not `Equal(a, b, eps)` — the
   arithmetic is identical, so exact equality must hold).
2. **`w_fext = 0` makes the constraint blind to the wrench.** Solve with
   `w_fext = 0` and the wrench applied; solve again with `w_fext = 0` and a
   zero wrench. Assert `constraint_torques` are equal within `eps` — the
   constraint solution must not depend on a wrench it is told to ignore.
3. **`w_fext = 1` makes the constraint fight the wrench.** Same two solves with
   `w_fext = 1`; assert the `constraint_torques` **differ** by more than `eps`
   on at least one joint. This is the contrast that gives case 2 meaning.
4. **Per-direction weighting is independent.** With
   `w_fext = (1, 0, 1)` (compliant only in the second constrained direction),
   assert the result differs from both the all-ones and the all-zeros results.
5. **Size validation.** `setDriverWeights` with a vector of size `nc + 1`
   returns `E_SIZE_MISMATCH`, and does not change the previously set weights
   (verify by re-solving and comparing to the pre-call result).

**Verify**: `cd orocos_kdl/build && make -j$(nproc) && ./tests/solvertest`
→ `OK (N tests)`, N one greater than before.

### Step 8: Full regression

**Verify**: `cd orocos_kdl/build && make check` → exit 0, every test passes,
including `VereshchaginTest` and `FdAndVereshchaginSolversConsistencyTest`
with unchanged expected values.

## Test plan

- New test `VereshchaginDriverWeightingTest` in
  `orocos_kdl/tests/solvertest.cpp`, five cases as listed in Step 7.
- Structural pattern: `SolverTest::VereshchaginTest()` at
  `orocos_kdl/tests/solvertest.cpp:818`.
- The single most important assertion is case 1's **exact** equality — it is
  what makes this change safe to merge into a stable release line.
- Verification: `cd orocos_kdl/build && make check` → all pass.

## Done criteria

ALL must hold:

- [ ] `cd orocos_kdl/build && make check` exits 0
- [ ] `VereshchaginTest` passes with its original hard-coded values, unedited
      (`git diff orocos_kdl/tests/solvertest.cpp | grep '669693'` → no output)
- [ ] `VereshchaginDriverWeightingTest` exists and all five cases pass
- [ ] `git diff orocos_kdl/src/chainhdsolver_vereshchagin.cpp` shows **no**
      change inside `final_upwards_sweep`
- [ ] `grep -n "using namespace Eigen" orocos_kdl/src/` → no output
- [ ] `CartToJnt`'s signature is unchanged
      (`git diff orocos_kdl/src/chainhdsolver_vereshchagin.hpp | grep 'int CartToJnt'` → no output)
- [ ] `git status --porcelain` lists only the four in-scope files
- [ ] `plans/README.md` status row updated

## STOP conditions

Stop and report back (do not improvise) if:

- Case 1 of Step 7 fails exact equality. That means the channel accumulators
  are not a clean decomposition and the "defaults are a no-op" guarantee is
  broken. Report which quantity diverges and by how much — do **not** relax
  the assertion to `eps`.
- `VereshchaginTest`'s expected values need editing to pass. They must not.
- You need to modify `final_upwards_sweep()`.
- Case 3 shows no difference — that would mean the solver is not compensating
  the wrench in the first place, contradicting this plan's premise.
- The `child.C` term appears in any `G_fext` / `G_ff` expression. It belongs
  to the unweighted nature part only; including it double-counts.

## Maintenance notes

- **Follow-ups deliberately excluded from this PR**, in likely order of value:
  (a) PyKDL binding for `setDriverWeights`; (b) the same channel split applied
  to `constraint_torques`, so a caller can ask for "the torque the wrench
  contributes" without a second solver instance; (c) weighting `beta` itself
  (the `w_ee` half of Eq. 3.42), which is currently achievable by the caller
  scaling `beta` before the call.
- Whether a weighted driver's joint torque must additionally be *commanded* is
  a caller-side question this API deliberately does not answer. A measured
  wrench (from an F/T sensor) physically exists and must not be commanded; a
  virtual wrench synthesised by a controller does not exist and must be. Both
  are legitimate; see plan 004 for how the downstream consumer handles it.
- A reviewer should check: the `(1 - w)` correction terms are exactly zero for
  default weights (not merely small); `G_fext` / `G_ff` never pick up `child.C`
  or the `v x (H v)` bias; and adding members to `segment_info` is acceptable
  for the target release line, since it changes the class layout.
- Interaction with plan 002: both edit `downwards_sweep` around lines 195–251.
  Land 002 first and rebase this on top.
