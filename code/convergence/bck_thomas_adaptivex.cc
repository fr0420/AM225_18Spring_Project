#include <fstream> // Header for file output
#include <vector> // Header for vectors
#include <cmath> // Header for sqrt(), fabs()
#include "omp.h"


/** Solves the fully irreversible CV process using backward Euler method with adaptive spatial grid with Thomas algorithm. */

int main() {
	printf("#Δx #rel_error #time_consumption\n");
	#pragma omp parallel for
    for(int power = 0;power<17;power++){
		// Specify simulation parameters
		double theta_i = 20.0;
		double theta_v = -20.0;
		double sigma = 100.0;
		//double deltaX = 2e-5;
		double deltaTheta = 0.02;
		// Specify the analytic peak flux
        double true_value = -0.446*sqrt(sigma);


		// Calculate other parameters
		double deltaT = deltaTheta / sigma;
		double maxT = 2 * fabs(theta_v - theta_i) / sigma;
		double maxX = 6 * sqrt(maxT);
		int m = (int)( maxT / deltaT ); // number of timesteps
		//printf("maxX=%g, deltaX=%g, n=%d, (n-1)*deltaX=%g\n", maxX, deltaX, n, (n-1)*deltaX);

		//unequally spaced grid
		double omega = 1.001;
		double h = 2e-5*pow(2,power);
		std::vector<double> X;
		X.push_back(0.0);
		while(X.back() < maxX) // where maxX = 6*sqrt(maxT)
		{
			X.push_back(X.back() + h);
			h *= omega;
		}
		int n = X.size(); //number of spacesteps
		//printf("maxX=%g, deltaX=%g, n=%d, (n-1)*deltaX=%g\n", maxX, deltaX, n, (n-1)*deltaX);

		// Calculate Thomas coefficients
		std::vector<double> alpha(n, 0);
		std::vector<double> beta(n, 0);
		std::vector<double> gamma(n, 0);
		for(int i=1; i<n-1; i++)
		{
			double delX_m = X[i] - X[i-1];
			double delX_p = X[i+1] - X[i];
			alpha[i] = -(2* deltaT) / ( delX_m * (delX_m + delX_p) );
			gamma[i] = -(2* deltaT) / ( delX_p * (delX_m + delX_p) );
			beta[i] = 1 - alpha[i] - gamma[i];
		}

        // Save voltage and current in an array
        double* data_array = new double[2*m];

		// Create containers
		std::vector<double> g_mod(n-1, 0);
		std::vector<double> C(n, 1.0); // concentration profile

		// Modify gamma coefficients
		// g_mod[0] = 0; // boundary condition
		for(int i=1; i<n-1; i++) {
			g_mod[i] = gamma[i] / (beta[i] - g_mod[i-1] * alpha[i]);
		}

		// Open file to output CV
		std::ofstream CV("CV_Output_changing_X.txt");

		// BEGIN SIMULATION
		double t0 = omp_get_wtime();
		double Theta = theta_i;
		for(int k=0; k<m; k++) {
			if ( k < m / 2 ) { 
				Theta -= deltaTheta; 
			} else { 
				Theta += deltaTheta; 
			}
		
			// Forward sweep - create modified deltas
			C[0] = 1.0 / (1.0 + exp(-Theta));
			for(int i=1; i<n-1; i++) {
				C[i] = ( C[i] - C[i-1]*alpha[i] ) / ( beta[i] - g_mod[i-1] * alpha[i] );
			}
			
			// Back Substitution
			// C[n-1] = 1.0;
			for(int i=n-2; i>=1; i--) {
				C[i] = C[i] - g_mod[i] * C[i+1];
			}
			
			// Output current
			// double flux = -(-C[2] + 4*C[1] -3*C[0]) / (2*(X[1]-X[0]));
			// printf("%g %g %g\n", omega, h, X[1]-X[0]);
			double flux = -( -C[2] + (omega+1.)*(omega+1.)*C[1] - omega*(omega+2.)*C[0] ) / ( (X[1]-X[0])*omega*(omega+1.) );
			CV << Theta << "\t" << flux << "\n";
	        data_array[k] = Theta;
            data_array[m+k]= flux;
		}
        //Calculate error
        //store the numerical peak flux
        double peak_flux = 0;

        for(int j = 0;j<m;j++){
            peak_flux = data_array[m+j]<peak_flux?data_array[m+j]:peak_flux;
        }
        double error = (peak_flux-true_value)/true_value;
    	double t1 = omp_get_wtime();
        double time_elapsed = t1-t0;
        //output dX, error and consumed time
        printf("%g %g %g\n",X[1]-X[0],fabs(error),time_elapsed);
	}
// END SIMULATION
}