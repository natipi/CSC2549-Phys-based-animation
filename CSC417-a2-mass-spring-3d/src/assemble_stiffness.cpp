#include <assemble_stiffness.h>

void assemble_stiffness(Eigen::SparseMatrixd &K, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXd> qdot, 
                     Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> E, Eigen::Ref<const Eigen::VectorXd> l0, 
                     double k) { 
        const rows = E.rows();
        int p0, p1;
        Eigen::MatrixXd Htemp(6,6); 
        f *= 0;

        typedef Eigen::Triplet<double> T;
        std::vector<T> tripletList;
        tripletList.reserve(estimation_of_entries);

        for(int i = 0; i < rows; i++) {
            p0 = E(i,0);
            p1 = E(i,1);

            // compute Hessian for single spring
            dV_spring_particle_dq(Htemp, q.segment(3*p0, 3*p0 + 2), q.segment(3*p1, 3*p1 + 2, l0(i), k(i)));

            // assemble. This ref https://cgvr.cs.uni-bremen.de/teaching/vr_1213/folien/10%20-%20mass%20spring%20systems-2.pdf helped clarify dimensions
            for(int k = 0; k < 6; k++) {
                for(int m = 0; m < 6; m++) {
                    tripletList.push_back(T((k < 3 ? 3 * p0 + k : 3 * p1 + k - 3),(m < 3 ? 3 * p0 + m : 3 * p1 + m - 3),Htemp(k,m)));
                }
            }
        }
        K.setFromTriplets(tripletList.begin(), tripletList.end());
    };