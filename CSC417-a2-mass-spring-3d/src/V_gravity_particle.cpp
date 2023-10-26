#include <V_gravity_particle.h>

// compute gravitational potential energy of a single particle
void V_gravity_particle(double &V, Eigen::Ref<const Eigen::Vector3d> q,  double mass, Eigen::Ref<const Eigen::Vector3d> g) {
    // V = mass * g * q(2); // get z coordinate for height
    V = mass * g.dot(q);
}