#include <d2V_spring_particle_particle_dq2.h>

// get the Hessian of V, -partial^2 V(q)/partial q^2
void d2V_spring_particle_particle_dq2(Eigen::Ref<Eigen::Matrix66d> H, Eigen::Ref<const Eigen::Vector3d> q0,  Eigen::Ref<const Eigen::Vector3d> q1, double l0, double stiffness) {
	float l = (qd - q0).norm();
    Eigen::MatrixXd Ipm = Matrix<double, 6,6>::Identity();
    Eigen::MatrixXd q(6, 1);

    for(int i = 0; i < 6; i++) {
        if (i < 3) {
            MatrixXd(i,0) = q0(i);
        }
        else {
            MatrixXd(i, 0) = -1 * q1(i-3);
            Ipm(i,i) = -1.0;
        }
    }

    H = stiffness * (l - l0)/l * Ipm - stiffness * q * q.transpose() / pow(l, 2);
}