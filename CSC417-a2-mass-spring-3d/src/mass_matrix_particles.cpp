#include <mass_matrix_particles.h>
#include <iostream>
#include <debug.h>

void mass_matrix_particles(Eigen::SparseMatrixd &M, Eigen::Ref<const Eigen::VectorXd> q, double mass) {
    const int len = q.size();
    M = Eigen::SparseMatrix<double> (len, len);
    
    for (int i = 0; i < len; i++){
        M.insert(i,i) = mass;
    }
}
