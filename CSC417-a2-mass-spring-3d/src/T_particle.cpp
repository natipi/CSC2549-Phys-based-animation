#include <T_particle.h>

// compute kinetic energy of a single particle
void T_particle(double &T, Eigen::Ref<const Eigen::VectorXd> qdot, double mass) {
    T = mass * qdot.norm();
}