#include <assemble_stiffness.h>
#include <iostream>

void assemble_stiffness(Eigen::SparseMatrixd &K, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXd> qdot, 
                     Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> E, Eigen::Ref<const Eigen::VectorXd> l0, 
                     double k) { 
        const int rows = E.rows();
        const int qlen = q.size();
        K = Eigen::SparseMatrixd (qlen, qlen);
        int p0, p1;
        Eigen::Matrix66d Htemp; 
        // K *= 0;

        typedef Eigen::Triplet<double> T;
        std::vector<T> tripletList;

        // estimated number of triples
        tripletList.reserve(36 * rows);

        for(int i = 0; i < rows; i++) {
            p0 = E(i,0);
            p1 = E(i,1);

            // compute Hessian for single spring
            d2V_spring_particle_particle_dq2(Htemp, q.segment(3*p0, 3), q.segment(3*p1, 3), l0(i), k);
            
            // assemble. This ref https://cgvr.cs.uni-bremen.de/teaching/vr_1213/folien/10%20-%20mass%20spring%20systems-2.pdf helped clarify dimensions
            for(int k = 0; k < 6; k++) {
                for(int m = 0; m < 6; m++) {
                    tripletList.push_back(T((k < 3 ? 3 * p0 + k : 3 * p1 + k - 3),(m < 3 ? 3 * p0 + m : 3 * p1 + m - 3),Htemp(k,m)));
                }
            }
        }
        K.setFromTriplets(tripletList.begin(), tripletList.end());
    };