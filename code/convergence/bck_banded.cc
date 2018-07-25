#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <vector>
#include <cstring>
#include <complex>   // std::complex
#include "omp.h"

typedef std::complex<double> dcomplex;

extern "C" {
    int zgbsv_(const int *n, const int *kl, const int *ku, const int *nrhs,
               dcomplex *ab, const int *ldab, int *ipiv, dcomplex *b,
               const int *ldb, int *info);
}


/** Solves the fully irreversible CV process using backward Euler method with banded solver. */

int main(){

	    // Specify simulation parameters
        double theta_i = 20.0;
        double theta_v = -20.0;
        double sigma = 100;
        // Store analytical result
        double true_value = -0.446*sqrt(sigma);
        
        double deltaTheta = 0.02;
        // Calculate other parameters
        double deltaT = deltaTheta / sigma;
        double maxT = 2 * fabs(theta_v - theta_i) / sigma;
        double maxX = 6*sqrt(maxT);
    #pragma omp parallel for
    for(int power = 0;power<17;power++){

        // DeltaX is a function of power
    	double deltaX = 2e-5*pow(2,power);
    	double deltaTheta = 0.02;
    	// Calculate other parameters
    	int n_space = (int)( maxX / deltaX ); // number of spacesteps
    	int m_time = (int)( maxT / deltaT ); // number of timesteps
    	// Calculate coefficients for banded solver
    	double lambda = deltaT / (deltaX*deltaX);
    	double alpha = -lambda;
    	double beta = 2.0*lambda + 1.0;
    	double gamma = -lambda;
        // Save voltage and current in an array
        double* data_array = new double[2*m_time];

//////////////////////////////////////////////////////////////////////////////////////////

    	// Use banded solver to solve the concentration profile

        // Size of the matrix
        const int n = n_space;
        // Number of lower and upper diagonals
        const int kl = 1, ku =1;
        // Number of source terms
        const int nrhs = 1;
        // Number of rows in banded storage
        const int ldab = 2*kl+ku+1;
        // Matrix in banded storage
        dcomplex* ab = new dcomplex[ldab*n];
        // Reuse ab
        dcomplex* ab_temp = new dcomplex[ldab*n];
        // Source term
        dcomplex* b = new dcomplex[n];
        // Pivot storage
        int* ipiv   = new int[n];
        // Error flag
        int info;
        // Set the dynamic potential 
    	double Theta = theta_i;

        // Matrix in row-major banded storage
        // Note that it is easier to construct the matrix on a per-diagonal basis
        // if we store things in row-major order
        std::vector<dcomplex> diags(ldab*n, 0.);

        // The first kl rows of banded storage are for internal use by LAPACK
        std::vector<dcomplex>::iterator start = diags.begin()+n*kl;

        // Construct the diagonals
        // Refer to LAPACK's banded storage format for details
        std::fill_n(start+n*kl,         n, dcomplex(beta,0));
        std::fill_n(start+n*(kl-1)+1, n-1, dcomplex(alpha,0));
        std::fill_n(start+n*(kl+1),   n-1, dcomplex(gamma,0));

        //Set the boundary conditions
        diags[n*2] = 1;
        diags[n+1] = 0;

        diags[n*2+n-1] = 1;
        diags[(ldab-1)*n+n-2] = 0;

        //Check the matrix construction
        /*
        for (int i=0; i<ldab; i++){
        	for(int j=0;j<n;j++){
        		printf("%g ",diags[n*i+j].real());
        	}
        	printf("\n");
        }*/

        // Transpose to column-major order for LAPACK
        for (int i=0; i<ldab; ++i)
            for (int j=0; j<n; ++j)
                ab[ldab*j+i] = diags[n*i+j];


        memcpy(ab_temp,ab,ldab*n*sizeof(dcomplex));
        for (int i=0; i<n; ++i) b[i] = dcomplex(1,0);
        //Boundary condition
        b[0] = dcomplex(1.0 / (1.0 + exp(-Theta)),0);
    	b[n-1] = dcomplex(1.0,0);
/////////////////////////////////////////////////////////////////////////////////
        //Perform simulation and count elapsed time
        double t0 = omp_get_wtime();
        for (int i=0 ; i<m_time;i++){

    		if(i<m_time/2) Theta -= deltaTheta; 
    		else Theta += deltaTheta; 
    	    b[0] = dcomplex(1.0 / (1.0 + exp(-Theta)),0);
    		b[n-1] = dcomplex(1.0,0);

    	    // Call LAPACK
    	    zgbsv_(&n, &kl, &ku, &nrhs, ab, &ldab, ipiv, b, &n, &info);

    		double flux = -(-b[2].real()+4*b[1].real()-3*b[0].real())/(2*deltaX);
            data_array[i] = Theta;
            data_array[m_time+i]= flux;
    
    		//reset ab (because zgbsv changes the composition of ab)
        	memcpy(ab,ab_temp,ldab*n*sizeof(dcomplex));

    	}
        //Calculate error
        //store the numerical peak flux
        double peak_flux = 0;

        for(int j = 0;j<m_time;j++){
            peak_flux = data_array[m_time+j]<peak_flux?data_array[m_time+j]:peak_flux;
        }

        double error = (peak_flux-true_value)/true_value;


    	double t1 = omp_get_wtime();
        double time_elapsed = t1-t0;
        //output dX, error and consumed time
        printf("%g %g %g\n",deltaX,fabs(error),time_elapsed);

        delete [] ab;
        delete [] ab_temp;
        delete [] b;
        delete [] ipiv;
        delete [] data_array;

    }
}