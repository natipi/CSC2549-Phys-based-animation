#include <fixed_point_constraints.h>
#include <algorithm>
#include <iostream>
#include <debug.h>

void fixed_point_constraints(Eigen::SparseMatrixd &P, unsigned int q_size, const std::vector<unsigned int> indices) {
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    tripletList.reserve(q_size);

    for (int i = 0; i < q_size; i++) {
        if(std::find(indices.begin(), indices.end(), i) != indices.end()) {
            // if i is in indices
            tripletList.push_back(T(i, i, 0));
        }
        else {
            tripletList.push_back(T(i,i, 1));
        }
    }

    P.setFromTriplets(tripletList.begin(), tripletList.end());

    debug("P nonzero entries: ");
    debug(P.nonZeros());
}