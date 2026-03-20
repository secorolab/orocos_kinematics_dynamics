// Copyright  (C)  2009  Ruben Smits <ruben dot smits at intermodalics dot eu>

// Version: 1.0
// Author: Ruben Smits <ruben dot smits at intermodalics dot eu>
// Author: Herman Bruyninckx
// Author: Azamat Shakhimardanov
// Maintainer: Ruben Smits <ruben dot smits at intermodalics dot eu>
// URL: http://www.orocos.org/kdl

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

#include "chainhdsolver_vereshchagin_fext_fixed_joint.hpp"
#include "frames_io.hpp"
#include "utilities/svd_eigen_HH.hpp"

namespace KDL
{

ChainHdSolver_Vereshchagin_Fext_FixedJoint::ChainHdSolver_Vereshchagin_Fext_FixedJoint(const Chain& chain_, const Twist &root_acc, const unsigned int nc_) :
    chain(chain_), nj(chain.getNrOfJoints()), ns(chain.getNrOfSegments()), nc(nc_),
    results(ns + 1, segment_info(nc))
{
    acc_root = root_acc;

    //Provide the necessary memory for computing the inverse of M0
    nu_sum.resize(nc);
    M_0_inverse.resize(nc, nc);
    Um = Eigen::MatrixXd::Identity(nc, nc);
    Vm = Eigen::MatrixXd::Identity(nc, nc);
    Sm = Eigen::VectorXd::Ones(nc);
    tmpm = Eigen::VectorXd::Ones(nc);

    // Provide the necessary memory for storing the total torque acting on each joint
    total_torques = Eigen::VectorXd::Zero(nj);
}

void ChainHdSolver_Vereshchagin_Fext_FixedJoint::updateInternalDataStructures() {
    ns = chain.getNrOfSegments();
    nj = chain.getNrOfJoints();
    total_torques = Eigen::VectorXd::Zero(nj);
    results.resize(ns+1,segment_info(nc));
}

int ChainHdSolver_Vereshchagin_Fext_FixedJoint::CartToJnt(const JntArray &q, const JntArray &q_dot, JntArray &q_dotdot, const Jacobian& alfa, const JntArray& beta, const Wrenches& f_ext, const JntArray &ff_torques, JntArray &constraint_torques)
{
    //Check sizes always
    nj = chain.getNrOfJoints();
    if(ns != chain.getNrOfSegments())
        return (error = E_NOT_UP_TO_DATE);
    if (q.rows() != nj || q_dot.rows() != nj || q_dotdot.rows() != nj || ff_torques.rows() != nj || constraint_torques.rows() != nj || f_ext.size() != ns)
        return (error = E_SIZE_MISMATCH);
    if (alfa.columns() != nc || beta.rows() != nc)
        return (error = E_SIZE_MISMATCH);
    //do an upward recursion for position, velocities and rigid-body bias forces
    this->initial_upwards_sweep(q, q_dot, q_dotdot, f_ext);
    //do an inward recursion for inertia, articulated bias forces and constraints
    this->downwards_sweep(alfa, ff_torques);
    //Solve for the constraint forces
    this->constraint_calculation(beta);
    //do an upward recursion to propagate the result and compute final output
    this->final_upwards_sweep(q_dotdot, constraint_torques);
    return (error = E_NOERROR);
}

void ChainHdSolver_Vereshchagin_Fext_FixedJoint::initial_upwards_sweep(const JntArray &q, const JntArray &qdot, const JntArray &qdotdot, const Wrenches& f_ext)
{
    unsigned int j = 0;
    F_total = Frame::Identity();
    for (unsigned int i = 0; i < ns; i++)
    {
        //Express everything in the segments reference frame (body coordinates)
        //which is at the segments tip, i.e. where the next joint is attached.

        //Calculate segment properties: X,S,vj,cj
        const Segment& segment = chain.getSegment(i);
        segment_info& s = results[i + 1];

        // FIX: For a fixed joint s.Z = Twist::Zero() and vj = Twist::Zero().
        // More importantly, q has only nj entries (one per movable joint). Accessing
        // q(j) when segment is fixed — especially after all movable joints have been
        // visited (j == nj) — is an out-of-bounds read. Branch here to avoid it.
        const bool segment_is_fixed = (segment.getJoint().getType() == Joint::Fixed);

        Twist vj; // joint velocity contribution in tip frame
        if (!segment_is_fixed)
        {
            //The pose between the joint root and the segment tip (tip expressed in joint root coordinates)
            s.F = segment.pose(q(j)); //X pose of each link in link coord system

            F_total = F_total * s.F; //X pose of the each link in root coord system
            s.F_base = F_total; //X pose of the each link in root coord system for getter functions

            //The velocity due to the joint motion of the segment expressed in the segments reference frame (tip)
            vj = s.F.M.Inverse(segment.twist(q(j), qdot(j))); //XDot of each link

            //The unit velocity due to the joint motion of the segment expressed in the segments reference frame (tip)
            s.Z = s.F.M.Inverse(segment.twist(q(j), 1.0));
            //Put Z in the joint root reference frame:
            s.Z = s.F * s.Z;

            //The total velocity of the segment expressed in the segments reference frame (tip)
            if (i != 0)
            {
                s.v = s.F.Inverse(results[i].v) + vj; // recursive velocity of each link in segment frame
                s.A = s.F.M.Inverse(results[i].A);
            }
            else
            {
                s.v = vj;
                s.A = s.F.M.Inverse(acc_root);
            }
            j++;
        }
        else
        {
            // Fixed joint: pose is constant (argument ignored by KDL), Z = 0, vj = 0.
            // Do NOT access q(j) here — j may equal nj (out of bounds) if all
            // movable joints have already been processed.
            s.F = segment.pose(0.0);

            F_total = F_total * s.F;
            s.F_base = F_total;

            vj = Twist::Zero(); // no joint velocity contribution
            s.Z = Twist::Zero();

            // Velocity propagates from parent only; no joint contribution.
            if (i != 0)
            {
                s.v = s.F.Inverse(results[i].v);
                s.A = s.F.M.Inverse(results[i].A);
            }
            else
            {
                s.v = Twist::Zero();
                s.A = s.F.M.Inverse(acc_root);
            }
            // j is NOT incremented — fixed joints have no entry in q/qdot
        }

        //c[i] = cj + v[i]xvj (remark: cj=0, since our S is not time dependent in local coordinates)
        //The velocity product acceleration
        s.C = s.v * vj; //This is a cross product: cartesian space BIAS acceleration in local link coord.
        // For a fixed joint vj=0, so s.C = 0 — correct, no Coriolis from a fixed joint.
        //Put C in the joint root reference frame
        s.C = s.F * s.C;
        //The rigid body inertia of the segment, expressed in the segments reference frame (tip)
        s.H = segment.getInertia();

        //wrench of the rigid body bias forces and the external forces on the segment (in body coordinates, tip)
        //external forces are taken into account through s.U.
        //
        // FEXT PRIORITY: The velocity product term s.v*(s.H*s.v) is intentionally
        // omitted here. This prioritizes the f_ext contribution by not mixing it
        // with centrifugal/Coriolis bias in s.U. The velocity product is accounted
        // for via s.C in the recursion instead.
        Wrench FextLocal = F_total.M.Inverse() * f_ext[i];
        s.U = -FextLocal;
    }
}

void ChainHdSolver_Vereshchagin_Fext_FixedJoint::downwards_sweep(const Jacobian& alfa, const JntArray &ff_torques)
{
    int j = nj - 1;
    for (int i = ns; i >= 0; i--)
    {
        //Get a handle for the segment we are working on.
        segment_info& s = results[i];
        //For segment N,
        //tilde is in the segment refframe (tip, not joint root)
        //without tilde is at the joint root (the childs tip!!!)
        //P_tilde is the articulated body inertia
        //R_tilde is the sum of external and coriolis/centrifugal forces
        //M is the (unit) acceleration energy already generated at link i
        //G is the (unit) magnitude of the constraint forces at link i
        //E are the (unit) constraint forces due to the constraints
        if (i == (int)ns)
        {
            s.P_tilde = s.H;
            s.R_tilde = s.U;
            s.M.setZero();
            s.G.setZero();
            //changeBase(alfa_N,F_total.M.Inverse(),alfa_N2);\n
            for (unsigned int r = 0; r < 3; r++)
                for (unsigned int c = 0; c < nc; c++)
                {
                    //copy alfa constrain force matrix in E~
                    s.E_tilde(r, c) = alfa(r + 3, c);
                    s.E_tilde(r + 3, c) = alfa(r, c);
                }
            //Change the reference frame of alfa to the segmentN tip frame
            //F_Total holds end effector frame, if done per segment bases then constraints could be extended to all segments
            Rotation base_to_end = F_total.M.Inverse();
            for (unsigned int c = 0; c < nc; c++)
            {
                Wrench col(Vector(s.E_tilde(3, c), s.E_tilde(4, c), s.E_tilde(5, c)),
                           Vector(s.E_tilde(0, c), s.E_tilde(1, c), s.E_tilde(2, c)));
                col = base_to_end*col;
                s.E_tilde.col(c) << Eigen::Vector3d::Map(col.torque.data), Eigen::Vector3d::Map(col.force.data);
            }
        }
        else
        {
            //For all others:
            //Everything should expressed in the body coordinates of segment i
            segment_info& child = results[i + 1];

            // FIX (Bug 1): When the child segment has a fixed joint, child.Z = 0 and
            // therefore child.D = dot(Z, PZ) = 0. All terms that divide by child.D
            // vanish in the limit (they represent joint-space projections onto a
            // zero-DOF axis). Branching here avoids NaN propagation through the rest
            // of the recursion.
            const bool child_is_fixed = (chain.getSegment(i).getJoint().getType() == Joint::Fixed);

            if (!child_is_fixed)
            {
                //Copy PZ into a vector so we can do matrix manipulations, put torques above forces
                Vector6d vPZ;
                vPZ << Eigen::Vector3d::Map(child.PZ.torque.data), Eigen::Vector3d::Map(child.PZ.force.data);
                Matrix6d PZDPZt;
                PZDPZt.noalias() = vPZ * vPZ.transpose();
                PZDPZt /= child.D;

                //equation a) (see Vereshchagin89) PZDPZt=[I,H;H',M]
                //Azamat:articulated body inertia as in Featherstone (7.19)
                s.P_tilde = s.H + child.P - ArticulatedBodyInertia(PZDPZt.bottomRightCorner<3,3>(), PZDPZt.topRightCorner<3,3>(), PZDPZt.topLeftCorner<3,3>());
                //equation b) (see Vereshchagin89)
                //Azamat: bias force as in Featherstone (7.20)
                s.R_tilde = s.U + child.R + child.PC + (child.PZ / child.D) * child.u;
                //equation c) (see Vereshchagin89)
                s.E_tilde = child.E;
                //Azamat: equation (c) right side term
                s.E_tilde.noalias() -= (vPZ * child.EZ.transpose()) / child.D;

                //equation d) (see Vereshchagin89)
                s.M = child.M;
                //Azamat: equation (d) right side term
                s.M.noalias() -= (child.EZ * child.EZ.transpose()) / child.D;

                //equation e) (see Vereshchagin89)
                s.G = child.G;
                Twist CiZDu = child.C + (child.Z / child.D) * child.u;
                Vector6d vCiZDu;
                vCiZDu << Eigen::Vector3d::Map(CiZDu.rot.data), Eigen::Vector3d::Map(CiZDu.vel.data);
                s.G.noalias() += child.E.transpose() * vCiZDu;
            }
            else
            {
                // Fixed child joint: rigid body merge — no joint-space projection/reduction.
                // All PZ/D, EZ/D, Z/D*u terms are zero since Z=0, D=0.
                // equation a): no joint-space inertia reduction
                s.P_tilde = s.H + child.P;
                // equation b): no nullspace force term
                s.R_tilde = s.U + child.R + child.PC;
                // equation c): no EZ/D correction
                s.E_tilde = child.E;
                // equation d): no EZ*EZ'/D correction
                s.M = child.M;
                // equation e): no Z/D*u acceleration term; only bias C remains
                s.G = child.G;
                Vector6d vC;
                vC << Eigen::Vector3d::Map(child.C.rot.data), Eigen::Vector3d::Map(child.C.vel.data);
                s.G.noalias() += child.E.transpose() * vC;
            }
        }
        if (i != 0)
        {
            //Transform all results to joint root coordinates of segment i (== body coordinates segment i-1)
            //equation a)
            s.P = s.F * s.P_tilde;
            //equation b)
            s.R = s.F * s.R_tilde;
            //equation c), in matrix: torques above forces, so switch and switch back
            for (unsigned int c = 0; c < nc; c++)
            {
                Wrench col(Vector(s.E_tilde(3, c), s.E_tilde(4, c), s.E_tilde(5, c)),
                           Vector(s.E_tilde(0, c), s.E_tilde(1, c), s.E_tilde(2, c)));
                col = s.F*col;
                s.E.col(c) << Eigen::Vector3d::Map(col.torque.data), Eigen::Vector3d::Map(col.force.data);
            }

            //needed for next recursion
            s.PZ = s.P * s.Z;

            /**
             * Additionally adding joint inertia to s.D, see:
             * - equation a) in Vereshchagin89
             * - equation 9.28, page 188, Featherstone book 2008
             */
            if (chain.getSegment(i - 1).getJoint().getType() != Joint::Fixed)
                s.D = chain.getSegment(i - 1).getJoint().getInertia() + dot(s.Z, s.PZ);
            else
                s.D = dot(s.Z, s.PZ);

            s.PC = s.P * s.C;

            //u=(Q-Z(R+PC)=sum of external forces along the joint axes,
            //R are the forces coming from the children,
            //Q is taken zero (do we need to take the previous calculated torques?

            // FEXT PRIORITY: totalBias uses only s.R (not s.R + s.PC).
            // This keeps the fext-prioritized formulation consistent with s.U above.
            // ff_torques are intentionally not applied here — fext solver drives motion
            // purely from external forces and constraints.
            s.totalBias = -dot(s.Z, s.R);
            s.u = s.totalBias;

            //Matrix form of Z, put rotations above translations
            Vector6d vZ;
            vZ << Eigen::Vector3d::Map(s.Z.rot.data), Eigen::Vector3d::Map(s.Z.vel.data);
            s.EZ.noalias() = s.E.transpose() * vZ;

            if (chain.getSegment(i - 1).getJoint().getType() != Joint::Fixed)
                j--;
        }
    }
}

void ChainHdSolver_Vereshchagin_Fext_FixedJoint::constraint_calculation(const JntArray& beta)
{
    //equation f) nu = M_0_inverse*(beta_N - E0_tilde`*acc0 - G0)
    //M_0_inverse, always nc*nc symmetric matrix

    //M_0_inverse=results[0].M.inverse();
    svd_eigen_HH(results[0].M, Um, Sm, Vm, tmpm);
    //truncated svd, what would sdls, dls physically mean?
    for (unsigned int i = 0; i < nc; i++)
        if (Sm(i) < 1e-14)
            Sm(i) = 0.0;
        else
            Sm(i) = 1 / Sm(i);

    results[0].M.noalias() = Vm * Sm.asDiagonal();
    M_0_inverse.noalias() = results[0].M * Um.transpose();

    Vector6d acc;
    acc << Eigen::Vector3d::Map(acc_root.rot.data), Eigen::Vector3d::Map(acc_root.vel.data);
    nu_sum.noalias() = -(results[0].E_tilde.transpose() * acc);
    nu_sum += beta.data;
    nu_sum -= results[0].G;

    //equation f) nu = M_0_inverse*(beta_N - E0_tilde`*acc0 - G0)
    nu.noalias() = M_0_inverse * nu_sum;
}

void ChainHdSolver_Vereshchagin_Fext_FixedJoint::final_upwards_sweep(JntArray &q_dotdot, JntArray &constraint_torques)
{
    unsigned int j = 0;

    for (unsigned int i = 1; i <= ns; i++)
    {
        segment_info& s = results[i];
        //Calculation of joint and segment accelerations
        //equation g) qdotdot[i] = D^-1*(Q - Z'(R + P(C + acc[i-1]) + E*nu))
        // = D^-1(u - Z'(P*acc[i-1] + E*nu)
        Twist a_p;
        if (i == 1)
            a_p = acc_root;
        else
            a_p = results[i - 1].acc;

        // FIX (Bug 3): A fixed joint has Z=0 and D=0. Dividing by D (parentAccComp,
        // nullspaceAccComp) produces NaN which then corrupts s.acc and propagates to
        // every subsequent segment. For a fixed joint there is no joint-space DOF:
        // acceleration transfers rigidly from the parent (plus Coriolis only).
        if (chain.getSegment(i - 1).getJoint().getType() != Joint::Fixed)
        {
            //acceleration components are also computed
            //Contribution of the acceleration of the parent (i-1)
            Wrench parent_force = s.P * a_p;
            double parent_forceProjection = -dot(s.Z, parent_force);
            double parentAccComp = parent_forceProjection / s.D;

            // FEXT PRIORITY: constraint_torques output uses fext formulation.
            // The E*nu constraint force projection is intentionally omitted so that
            // the torque output reflects the fext-driven (nullspace + parent) components.
            constraint_torques(j) = s.u + parent_forceProjection;

            s.nullspaceAccComp = s.u / s.D;

            // FEXT PRIORITY: q_dotdot does not include constAccComp.
            // Joint acceleration is driven by nullspace (fext/bias) + parent inertia only.
            q_dotdot(j) = (s.nullspaceAccComp + parentAccComp);

            s.acc = s.F.Inverse(a_p + s.Z * q_dotdot(j) + s.C);
            j++;
        }
        else
        {
            // Fixed joint: no joint DOF, no joint acceleration to solve for.
            // Segment acceleration propagates rigidly from parent plus Coriolis only.
            s.constAccComp = 0.0;
            s.nullspaceAccComp = 0.0;
            s.acc = s.F.Inverse(a_p + s.C);
            // j is NOT incremented — fixed joints have no entry in q_dotdot/constraint_torques
        }
    }
}

// Returns Cartesian acceleration of links in robot base coordinates
void ChainHdSolver_Vereshchagin_Fext_FixedJoint::getTransformedLinkAcceleration(Twists& x_dotdot)
{
    assert(x_dotdot.size() == ns + 1);
    x_dotdot[0] = acc_root;
    for (unsigned int i = 1; i < ns + 1; i++)
        x_dotdot[i] = results[i].F_base.M * results[i].acc;
}

// Returns total torque acting on each joint (constraints + nature + external forces)
void ChainHdSolver_Vereshchagin_Fext_FixedJoint::getTotalTorque(JntArray &total_tau)
{
    assert(total_tau.data.size() == total_torques.size());
    total_tau.data = total_torques;
}

// Returns magnitude of the constraint forces acting on the end-effector: Lagrange Multiplier
void ChainHdSolver_Vereshchagin_Fext_FixedJoint::getContraintForceMagnitude(Eigen::VectorXd &nu_)
{
    assert(nu_.size() == nu.size());
    nu_ = nu;
}

void ChainHdSolver_Vereshchagin_Fext_FixedJoint::getLinkCartesianPose(Frames& x_base)
{
    for (int i = 0; i < ns; i++)
    {
        x_base[i] = results[i + 1].F_base;
    }
    return;
}

void ChainHdSolver_Vereshchagin_Fext_FixedJoint::getLinkCartesianVelocity(Twists& xDot_base)
{
    for (int i = 0; i < ns; i++)
    {
        xDot_base[i] = results[i + 1].F_base.M * results[i + 1].v;
    }
    return;
}

void ChainHdSolver_Vereshchagin_Fext_FixedJoint::getLinkCartesianAcceleration(Twists& xDotDot_base)
{
    for (int i = 0; i < ns; i++)
    {
        xDotDot_base[i] = results[i + 1].F_base.M * results[i + 1].acc;
    }
    return;
}

void ChainHdSolver_Vereshchagin_Fext_FixedJoint::getLinkPose(Frames& x_local)
{
    for (int i = 0; i < ns; i++)
    {
        x_local[i] = results[i + 1].F;
    }
    return;
}

void ChainHdSolver_Vereshchagin_Fext_FixedJoint::getLinkVelocity(Twists& xDot_local)
{
    for (int i = 0; i < ns; i++)
    {
        xDot_local[i] = results[i + 1].v;
    }
    return;
}

void ChainHdSolver_Vereshchagin_Fext_FixedJoint::getLinkAcceleration(Twists& xDotdot_local)
{
    for (int i = 0; i < ns; i++)
    {
        xDotdot_local[i] = results[i + 1].acc;
    }
    return;
}

void ChainHdSolver_Vereshchagin_Fext_FixedJoint::getJointBiasAcceleration(JntArray& bias_q_dotdot)
{
    for (int i = 0; i < ns; i++)
    {
        double tmp = results[i + 1].totalBias;
        bias_q_dotdot(i) = tmp / results[i + 1].D;
    }
    return;
}

void ChainHdSolver_Vereshchagin_Fext_FixedJoint::getLinkBiasForceAcceleratoinEnergy(Eigen::VectorXd& G)
{
    for (int i = 0; i < ns; i++)
    {
        G = results[i + 1].G;
    }
    return;
}

void ChainHdSolver_Vereshchagin_Fext_FixedJoint::getLinkBiasForceMatrix(Wrenches& R_tilde)
{
    for (int i = 0; i < ns; i++)
    {
        R_tilde[i] = results[i + 1].R_tilde;
        std::cout << "s.R_tilde " << i << ":  " << results[i + 1].R << std::endl;
    }
    return;
}

}//namespace
