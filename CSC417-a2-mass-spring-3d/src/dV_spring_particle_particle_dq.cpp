#include <dV_spring_particle_particle_dq.h>

// Compute the forces exerted by a single spring on its end-points.
void dV_spring_particle_particle_dq(Eigen::Ref<Eigen::Vector6d> f, Eigen::Ref<const Eigen::Vector3d> q0,  Eigen::Ref<const Eigen::Vector3d> q1, double l0, double stiffness) {
    Eigen::MatrixXd C(6, 3);
    C << 1, 0, 0, 
         0, 1, 0, 
         0, 0, 1,
         -1, 0, 0,
         0, -1, 0,
         0, 0, -1;
    f = C * (stiffness * ((q1 - q0).norm() - l0) * (q1 - q0) / (q1 - q0).norm());
}