#include <mass_matrix_particles.h>

// TO DO: WHY DO I NEED Q???
void mass_matrix_particles(Eigen::SparseMatrixd &M, Eigen::Ref<const Eigen::VectorXd> q, double mass) {
    M.diagonal().array() += mass;
}
