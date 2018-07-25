#include "butler_volmer.hh"

#include <sstream>
#include <cmath>
#include "omp.h"


/** Main experiments done here. */
int main() {

	const double theta_i = 20.0;
	const double theta_v = -20.0;
	const double sigma = 100.0;
	const double deltaTheta = 0.02;
	const double k_not = 50.;
	// const double trans_coeff = 0.5;
	const double omega = 1.001;
	const double h = 0.00016;

	double k_list[] = {900., 50., 2., 0.15625};

#pragma omp parallel for
	// for(int i=0; i<4; i++) {
	// 	double k_not = k_list[i];
	// for(int p=-10; p<8; p++){
	// 	double k_not = 5*pow(2,p); 
	for(int l=-4; l<5; l++) {
		double trans_coeff = 0.5+l*0.1;
		std::string s1 = "./snaps_k=50/alpha=", s2 = ".out";
		std::stringstream ss;
  		ss << s1 << trans_coeff << s2;
  		const std::string str = ss.str();
  		const char *fn = str.c_str();

	    butler_volmer bv(theta_i,theta_v,sigma,deltaTheta,k_not,trans_coeff);
	    bv.init_grid(omega,h);
	    bv.solve(fn,4);
	}
}
