// Probe: what does each ACHD variant actually return, and is the two-pass sum a clean split?
// Diagnostic only, not part of the library. Evidence for plans/004.
//
//   g++ -O0 -g -std=c++17 plans/reference/achd_probe.cpp \
//       -I <ws>/install/orocos_kdl/include -I /usr/include/eigen3 \
//       -L <ws>/install/orocos_kdl/lib -lorocos-kdl -Wl,-rpath,<ws>/install/orocos_kdl/lib -o /tmp/achd_probe
//
// Shows: the _fext pass returns ~gravity torque with a ZERO wrench, and still
// returns torque from joint velocity alone. It is not a clean f_ext channel.
// Requires the forked _fixed_joint solvers, so it stops compiling once
// plans/004 deletes them -- that is intended, it is a pre-migration probe.
#include <kdl/chain.hpp>
#include <kdl/chainhdsolver_vereshchagin_fixed_joint.hpp>
#include <kdl/chainhdsolver_vereshchagin_fext_fixed_joint.hpp>
#include <cstdio>

using namespace KDL;

static Chain make_arm(int n, double len, double mass)
{
    Chain c;
    for (int i = 0; i < n; i++)
    {
        RigidBodyInertia I(mass, Vector(0, 0, len / 2),
                           RotationalInertia(mass * len * len / 12.0, mass * len * len / 12.0, 1e-4));
        c.addSegment(Segment(Joint(Joint::RotY), Frame(Vector(0, 0, len)), I));
    }
    return c;
}

static void dump(const char* tag, const JntArray& t)
{
    printf("%-28s", tag);
    for (unsigned i = 0; i < t.rows(); i++) printf("% 10.4f", t(i));
    printf("\n");
}

int main()
{
    const int N = 3;
    Chain chain = make_arm(N, 0.4, 2.0);
    Twist root_acc(Vector(0, 0, 9.81), Vector::Zero()); // KDL ACHD convention, as generated code uses

    const unsigned nc = 6;
    ChainHdSolver_Vereshchagin_Fixed_Joint acc(chain, root_acc, nc);
    ChainHdSolver_Vereshchagin_Fext_FixedJoint fext(chain, root_acc, nc);

    JntArray q(N), qd(N), qdd(N), ff(N), tau_acc(N), tau_fext(N);
    q(0) = 0.3; q(1) = -0.6; q(2) = 0.4;   // non-trivial pose
    qd(0) = 0.0; qd(1) = 0.0; qd(2) = 0.0; // start at rest: isolates gravity

    Jacobian alpha(nc), alpha_zero(nc);
    JntArray beta(nc), beta_zero(nc), ff_zero(N);
    SetToZero(alpha_zero); SetToZero(beta_zero); SetToZero(ff_zero); SetToZero(ff);
    alpha.data.setIdentity();          // full 6-DoF constraint
    SetToZero(beta);                   // hold still: zero acceleration energy

    Wrenches f_ext(chain.getNrOfSegments());
    Wrenches f_ext_zero(chain.getNrOfSegments());

    printf("\n--- at rest (qd=0), NO external wrench -------------------------------\n");
    fext.CartToJnt(q, qd, qdd, alpha_zero, beta_zero, f_ext_zero, ff_zero, tau_fext);
    acc.CartToJnt(q, qd, qdd, alpha, beta, f_ext_zero, ff, tau_acc);
    dump("tau_fext (should be ~0)", tau_fext);
    dump("tau_acc  (gravity comp)", tau_acc);

    printf("\n--- at rest, 10 N downward wrench on the last segment ----------------\n");
    f_ext[chain.getNrOfSegments() - 1] = Wrench(Vector(0, 0, -10.0), Vector::Zero());
    fext.CartToJnt(q, qd, qdd, alpha_zero, beta_zero, f_ext, ff_zero, tau_fext);
    acc.CartToJnt(q, qd, qdd, alpha, beta, f_ext_zero, ff, tau_acc);
    dump("tau_fext (wrench pass)", tau_fext);
    dump("tau_acc  (accel pass)", tau_acc);
    JntArray sum(N);
    Add(tau_fext, tau_acc, sum);
    dump("sum (what is commanded)", sum);

    printf("\n--- reference: single solver, both drivers at once -------------------\n");
    JntArray tau_both(N);
    acc.CartToJnt(q, qd, qdd, alpha, beta, f_ext, ff, tau_both);
    dump("tau_acc(f_ext given)", tau_both);

    printf("\n--- moving (qd != 0), no wrench: Coriolis leak check -----------------\n");
    qd(0) = 1.0; qd(1) = -1.5; qd(2) = 2.0;
    fext.CartToJnt(q, qd, qdd, alpha_zero, beta_zero, f_ext_zero, ff_zero, tau_fext);
    dump("tau_fext (should be ~0)", tau_fext);

    // Band-aid candidate: give the wrench pass a zero root acceleration, so it
    // carries no gravity. Does the remaining output become the pure f_ext term?
    printf("\n--- fext pass built with root_acc = 0 --------------------------------\n");
    ChainHdSolver_Vereshchagin_Fext_FixedJoint fext0(chain, Twist::Zero(), nc);
    JntArray tau_f0_rest(N), tau_f0_move(N), tau_f0_wrench(N);
    SetToZero(qd);
    fext0.CartToJnt(q, qd, qdd, alpha_zero, beta_zero, f_ext_zero, ff_zero, tau_f0_rest);
    dump("at rest, no wrench", tau_f0_rest);
    fext0.CartToJnt(q, qd, qdd, alpha_zero, beta_zero, f_ext, ff_zero, tau_f0_wrench);
    dump("at rest, 10N wrench", tau_f0_wrench);
    qd(0) = 1.0; qd(1) = -1.5; qd(2) = 2.0;
    fext0.CartToJnt(q, qd, qdd, alpha_zero, beta_zero, f_ext_zero, ff_zero, tau_f0_move);
    dump("moving, no wrench", tau_f0_move);
    return 0;
}
