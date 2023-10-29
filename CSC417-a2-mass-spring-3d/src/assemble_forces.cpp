#include <assemble_forces.h>
#include <iostream>

// Iterate through each spring in the mesh, compute the per-spring forces and assemble the global force vector.
void assemble_forces(Eigen::VectorXd &f, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXd> qdot, 
                     Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> E, Eigen::Ref<const Eigen::VectorXd> l0, 
                     double mass, double k) { 
        const int rows = E.rows();
        int p0, p1;
        Eigen::Vector6d ftemp; 
        f = Eigen::VectorXd (q.size());
        f *= 0;
        for(int i = 0; i < rows; i++) {
            p0 = E(i,0);
            p1 = E(i,1);

            // compute 6d force vector (f, -f) force on start and end spring
            dV_spring_particle_particle_dq(ftemp, q.segment(3*p0, 3), q.segment(3*p1, 3), l0(i), k);

            // I have to += cause im adding the forces from every spring. But do I need to initialize f = 0 since it might have some value from before? 
            f.segment(3*p0, 3) += ftemp.segment(0,3);
            f.segment(3*p1, 3) += ftemp.segment(3,3);
        }
    };