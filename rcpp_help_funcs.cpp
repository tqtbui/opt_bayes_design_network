// [[Rcpp::depends(RcppEigen)]]

#include <RcppEigen.h>

// [[Rcpp::export]]
SEXP cpp_pow_crit(
	const Eigen::Map<Eigen::MatrixXd> m, Eigen::Map<Eigen::MatrixXd> m_der, 
	Eigen::Map<Eigen::MatrixXd> beta, 
	double netfunc, double netfunc_der, double N, double omega) {

	Eigen::MatrixXd fish = Eigen::MatrixXd::Zero(6, 6);
	fish.topLeftCorner(4, 4) = (1/omega) * m.adjoint() * m;
	fish(5, 5) = N / ( 2 * pow(omega, 2) );
	fish(4, 4) = (1/omega) * (beta.transpose() * m_der.transpose() * m_der * beta).value();
	fish.block<4,1>(0,4) = (1/omega) * m.adjoint() * m_der * beta;

	fish.triangularView<Eigen::Lower>() = fish.transpose();
	Eigen::MatrixXd hess = fish.inverse();

	Eigen::MatrixXd dtotal = Eigen::MatrixXd::Zero(5, 1);
	dtotal(1, 0) = 1;
	dtotal(2, 0) = netfunc; 
	dtotal(3, 0) = - netfunc;
	dtotal(4, 0) = (beta(2,0) - beta(3,0))*netfunc_der;

	double ret = (dtotal.adjoint() * hess.topLeftCorner(5, 5) * dtotal).value();
	return Rcpp::wrap(ret);
}

// [[Rcpp::export]]
SEXP cpp_locopt_crit(
	const Eigen::Map<Eigen::MatrixXd> D, Eigen::MappedSparseMatrix<double> W, 
	double rho, Eigen::Map<Eigen::MatrixXd> x) {

	Eigen::MatrixXd fish = x.adjoint() * D * x - rho * x.adjoint() * W * x;
	Eigen::MatrixXd hess = fish.inverse();

	return Rcpp::wrap(hess(1,1));
}

// [[Rcpp::export]]
SEXP cpp_basse_crit(
	const Eigen::MappedSparseMatrix<double> AA, 
	Eigen::Map<Eigen::VectorXd> d, double N,
	double mu, double gamma, double sigma, Eigen::Map<Eigen::VectorXd> z1) {

	double n1 = z1.sum();
	Eigen::VectorXd z0 = Eigen::VectorXd::Ones(N) - z1;
	double n0 = z0.sum();
	Eigen::VectorXd AAz1 = AA * z1;

	double bias = (1/n1) * z1.dot(d) - (1/n0) * z0.dot(d);
	double var1 = (1/n1) + (1/n0);
	double var2 = (1/pow(n1, 2)) * z1.dot(AAz1) + (1/pow(n0, 2)) * z0.dot(AA * z0) - (2/(n1*n0)) * z0.dot(AAz1);

	double ret = pow(mu, 2) * pow(bias, 2) + pow(gamma, 2) * var1 + pow(sigma, 2) * var2;

	return Rcpp::wrap(ret);
}