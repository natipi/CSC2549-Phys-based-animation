//Input:
//  q - generalized coordiantes for the mass-spring system
//  qdot - generalized velocity for the mass spring system
//  dt - the time step in seconds
//  mass - the mass
//  force(q, qdot) - a function that computes the force acting on the mass as a function. This takes q and qdot as parameters.
//Output:
//  q - set q to the updated generalized coordinate using Runge-Kutta time integration
//  qdot - set qdot to the updated generalized velocity using Runge-Kutta time integration

template<typename FORCE> 
inline void runge_kutta(Eigen::VectorXd &q, Eigen::VectorXd &qdot, double dt, double mass,  FORCE &force) {
	// store slopes kappa1, kappa2, each of which themselves has a q and qdot coordinate
	// recall A^-1 is 1/m in the qdot coord and 1 in the q coord
	Eigen::VectorXd f;
	Eigen::VectorXd qprime; 
	Eigen::VectorXd qdotprime; 

	// k1 = f(y(t0), t0)
	// f = ( 0 & -k/m \\ 1 & 0 ) so sends (v, q) to (-kq/m, v)
	force(f, q, qdot); // the content of f is now -kq, namely first coordinate of f(\vec y)
	Eigen::VectorXd k10 = f / mass;
	Eigen::VectorXd k11 = qdot; 

	// k2 = f(y(t0) + k1 dt/2, t0 + dt/2)
	qdotprime = qdot + dt * k10 / 2;
	qprime = q + dt * k11 / 2;
	force(f, qprime, qdotprime);
	Eigen::VectorXd k20 = f / mass;
	Eigen::VectorXd k21 = qdotprime;

	qprime = q + dt * k21 / 2;
	qdotprime = qdot + dt * k20 / 2;
	force(f, qprime, qdotprime);
	Eigen::VectorXd k30 = f / mass;
	Eigen::VectorXd k31 = qdotprime;

	qprime = q + dt * k31; 
	qdotprime = qdot + dt * k30;
	force(f, qprime, qdotprime);
	Eigen::VectorXd k40 = f / mass;
	Eigen::VectorXd k41 = qdotprime;

	qdot += dt * (k10 + 2 * k20 + 2 * k30 + k40) / 6;
	q += dt * (k11 + 2 * k21 + 2 * k31 + k41) / 6; 

}