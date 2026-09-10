// Probe: does `beta` mean the TRUE end-effector acceleration, or the gravity-offset
// pseudo-acceleration? Compares pristine upstream against the secoro fork's E_input term.
// Diagnostic only, not part of the library. Evidence for plans/001 and plans/004.
//
// Extract a renamed copy of the pristine upstream solver, then build:
//
//   rm -rf /tmp/upver && mkdir -p /tmp/upver
//   git show origin/master:orocos_kdl/src/chainhdsolver_vereshchagin.hpp > /tmp/upver/up.hpp
//   git show origin/master:orocos_kdl/src/chainhdsolver_vereshchagin.cpp > /tmp/upver/up.cpp
//   sed -i 's/ChainHdSolver_Vereshchagin/ChainHdSolver_Vereshchagin_Up/g; s/KDL_CHAINHDSOLVER_VERESHCHAGIN_HPP/KDL_UPVER_HPP/g' /tmp/upver/up.hpp /tmp/upver/up.cpp
//   sed -i 's|#include "chainhdsolver_vereshchagin.hpp"|#include "up.hpp"|' /tmp/upver/up.cpp
//   sed -i 's|#include "\(chainidsolver\|frames\|articulatedbodyinertia\).hpp"|#include <kdl/\1.hpp>|' /tmp/upver/up.hpp
//   sed -i 's|#include "frames_io.hpp"|#include <kdl/frames_io.hpp>|; s|#include "utilities/svd_eigen_HH.hpp"|#include <kdl/utilities/svd_eigen_HH.hpp>|' /tmp/upver/up.cpp
//
//   g++ -O0 -g -std=c++17 plans/reference/achd_beta_semantics.cpp /tmp/upver/up.cpp \
//       -I <ws>/install/orocos_kdl/include -I /usr/include/eigen3 -I /tmp/upver \
//       -L <ws>/install/orocos_kdl/lib -lorocos-kdl -Wl,-rpath,<ws>/install/orocos_kdl/lib -o /tmp/beta_sem
//
// NOTE: after plans/004 removes E_input, the two columns converge and the
// "fork" row stops differing -- that is the expected post-migration result.
#include <kdl/chain.hpp>
#include <kdl/chainhdsolver_vereshchagin.hpp>   // fork build (has E_input)
#include "/tmp/upver/up.hpp"                    // pristine upstream (no E_input)
#include <cstdio>
#include <cmath>
#include <algorithm>

using namespace KDL;

static Chain make_planar_arm(int n, double len, double mass)
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

static void report(const char* tag, const JntArray& tau, const JntArray& qdd, const Twist& ee)
{
    printf("  %-22s tau=[% 8.3f % 8.3f % 8.3f]  qdd=[% 8.3f % 8.3f % 8.3f]  ee_acc=(x % 7.3f, z % 7.3f)\n",
           tag, tau(0), tau(1), tau(2), qdd(0), qdd(1), qdd(2), ee.vel.x(), ee.vel.z());
}

int main()
{
    const int N = 3;
    Chain chain = make_planar_arm(N, 0.4, 2.0);
    const unsigned ns = chain.getNrOfSegments();
    Twist root_acc(Vector(0, 0, 9.81), Vector::Zero());

    const unsigned nc = 3;
    Jacobian alpha(nc);
    alpha.setColumn(0, Twist(Vector(1, 0, 0), Vector::Zero())); // linear x
    alpha.setColumn(1, Twist(Vector(0, 0, 1), Vector::Zero())); // linear z
    alpha.setColumn(2, Twist(Vector::Zero(), Vector(0, 1, 0))); // angular y

    JntArray q(N), qd(N), qdd(N), ff(N), tau(N);
    q(0) = 0.3; q(1) = -0.6; q(2) = 0.4;
    SetToZero(qd); SetToZero(ff);
    Wrenches f_ext(ns);
    std::vector<Twist> xdd(ns + 1);

    ChainHdSolver_Vereshchagin_Up up(chain, root_acc, nc);
    ChainHdSolver_Vereshchagin fork(chain, root_acc, nc);

    JntArray beta(nc);

    printf("\nTask: hold the end-effector still  =>  true EE acceleration must be 0\n");

    printf("\n[beta = 0]  (what a user writes, per the header doc `beta = alpha^T * Xdd_N`)\n");
    SetToZero(beta);
    up.CartToJnt(q, qd, qdd, alpha, beta, f_ext, ff, tau);
    up.getTransformedLinkAcceleration(xdd);
    report("upstream", tau, qdd, xdd[ns]);
    fork.CartToJnt(q, qd, qdd, alpha, beta, f_ext, ff, tau);
    fork.getTransformedLinkAcceleration(xdd);
    report("fork (E_input)", tau, qdd, xdd[ns]);

    printf("\n[beta = alpha^T * root_acc = (0, +9.81, 0)]  (caller-side shift, no solver change)\n");
    SetToZero(beta);
    beta(1) = 9.81; // linear z column
    up.CartToJnt(q, qd, qdd, alpha, beta, f_ext, ff, tau);
    up.getTransformedLinkAcceleration(xdd);
    report("upstream + shift", tau, qdd, xdd[ns]);

    // Equivalence check across a sweep of poses and a non-zero beta: is the fork's
    // E_input term identical to passing beta + alpha^T*root_acc to upstream?
    printf("\nEquivalence sweep (fork[beta] vs upstream[beta + alpha^T*root_acc]):\n");
    double worst = 0.0;
    JntArray tau_u(N), tau_f(N), qdd_u(N), qdd_f(N), b_f(nc), b_u(nc);
    for (int k = 0; k < 12; k++)
    {
        q(0) = -1.2 + 0.21 * k; q(1) = 0.9 - 0.17 * k; q(2) = 0.4 + 0.11 * k;
        qd(0) = 0.5 * k - 1.0; qd(1) = 1.3 - 0.2 * k; qd(2) = -0.7 + 0.3 * k;
        ff(0) = 1.5; ff(1) = -2.0; ff(2) = 0.5;
        f_ext[ns - 1] = Wrench(Vector(3.0, 0.0, -7.0), Vector(0.0, 1.0, 0.0));
        b_f(0) = 0.3; b_f(1) = -1.1; b_f(2) = 0.05;
        b_u(0) = b_f(0); b_u(1) = b_f(1) + 9.81; b_u(2) = b_f(2); // + alpha^T*root_acc
        fork.CartToJnt(q, qd, qdd_f, alpha, b_f, f_ext, ff, tau_f);
        up.CartToJnt(q, qd, qdd_u, alpha, b_u, f_ext, ff, tau_u);
        for (int i = 0; i < N; i++)
        {
            worst = std::max(worst, std::abs(tau_u(i) - tau_f(i)));
            worst = std::max(worst, std::abs(qdd_u(i) - qdd_f(i)));
        }
    }
    printf("  max |difference| in tau and qdd over 12 poses: %.3e\n", worst);

    printf("\nNote: ee_acc is what getTransformedLinkAcceleration() reports, which upstream's\n"
           "VereshchaginTest asserts equal to beta. Check which row satisfies that.\n\n");
    return 0;
}
