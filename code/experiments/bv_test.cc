#include "butler_volmer.hh"

/** A minimal example showing how to use the butler_volmer solver. */
int main() {

	const double theta_i = 20.0;
	const double theta_v = -20.0;
	const double sigma = 100.0;
	const double deltaTheta = 0.02;
	const double k_not = 640;
	const double trans_coeff = 0.5;
	const double omega = 1.001;
	const double h = 0.00016;

    butler_volmer bv(theta_i,theta_v,sigma,deltaTheta,k_not,trans_coeff);
    bv.init_grid(omega,h);
    bv.solve("test_snaps.out",4);
    bv.output_cv("test_cv.out");
    bv.print_peaks();

}
