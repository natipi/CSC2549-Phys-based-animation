#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <EigenTypes.h>
#include <iostream>
#include <debug.h>

//Input:
//  q - generalized coordiantes for the mass-spring system
//  qdot - generalized velocity for the mass spring system
//  dt - the time step in seconds
//  mass - the mass matrix
//  force(q, qdot) - a function that computes the force acting on the mass-spring system. This takes q and qdot as parameters.
//  stiffness(q, qdot) - a function that computes the stiffness (negative second derivative of the potential energy). This takes q and qdot as parameters.  
//  tmp_force - scratch space to collect forces
//  tmp_stiffness - scratch space to collect stiffness matrix
//Output:
//  q - set q to the updated generalized coordinate using linearly implicit time integration
//  qdot - set qdot to the updated generalized velocity using linearly implicit time integration
template<typename FORCE, typename STIFFNESS> 
inline void linearly_implicit_euler(Eigen::VectorXd &q, Eigen::VectorXd &qdot, double dt, 
                            const Eigen::SparseMatrixd &mass,  FORCE &force, STIFFNESS &stiffness, 
                            Eigen::VectorXd &tmp_force, Eigen::SparseMatrixd &tmp_stiffness) {
    stiffness(tmp_stiffness, q, qdot);
    force(tmp_force, q, qdot);

    debug2(tmp_force.norm());

    // Eigen::SparseLU<Eigen::SparseMatrix<double> > solver;
    
    Eigen::SparseMatrix<double> a = mass - dt * dt * tmp_stiffness;
    Eigen::SimplicialCholesky<Eigen::SparseMatrix<double> > solver(a);

    // debug2("dt ");debug2(dt);
    // debug2("mass nonzeroes");debug2(mass.nonZeros());
    // debug2("stiffness nonzeroes");debug2(tmp_stiffness.nonZeros());
    // int massrows = mass.rows();
    // std::cout << "mass diagonal: ";
    // for (int i = 0; i < massrows; i++) {
    //     std::cout << mass.coeff(i,i) << " ";
    // }
    // debug(" ");
    // a.makeCompressed();
    // std::cout << "a matrix: " << a << std::endl << std::endl << std::endl;
    
    // solver.compute(a);

    // if(solver.info()!=Eigen::Success) {
    //     if (solver.info()==Eigen::NumericalIssue)
    //         debug2("numerical issue");
    //     if (solver.info()==Eigen::InvalidInput)
    //         debug2("invalid input ");
    //     return;
    // }

    qdot = solver.solve(mass * qdot + dt * tmp_force);

    // if(solver.info()!=Eigen::Success) {
    //     debug2("solve failed");
    //     return;
    // }

    q += dt * qdot; 
}