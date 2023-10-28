#include <fixed_point_constraints.h>
#include <algorithm>

void fixed_point_constraints(Eigen::SparseMatrixd &P, unsigned int q_size, const std::vector<unsigned int> indices) {
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    tripletList.reserve(estimation_of_entries);

    int row = 0;
    for (int i = 0; i < q_size; i++) {
        if(std::find(indices.begin(), indices.end(), i) == indices.end()) {
            // if i is not in indices
            tripletList.push_back(T(row, i, 1));
            row++;
        }
    }

    P.setFromTriplets(tripletList.begin(), tripletList.end());
}