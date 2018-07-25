#include "butler_volmer.hh"

#include <cstdio>
#include <cstring>
#include <vector> // Header for vectors
#include <cmath> // Header for sqrt(), fabs()

/** The class constructor sets constants and dynamically allocate memory. 
 * \param[in] theta_i_ initial potential.
 * \param[in] theta_v_ vertex potential.
 * \param[in] sigma_ scan rate.
 * \param[in] dtheta_ size of potential step.
 * \param[in] k0_ standard heterogeneous rate.
 * \param[in] trans_ charge transfer coefficient. 
 */
butler_volmer::butler_volmer(double theta_i_,double theta_v_,double sigma_,double dtheta_,double k0_,double trans_) : 
    theta_i(theta_i_), theta_v(theta_v_), sigma(sigma_), dtheta(dtheta_), k0(k0_), trans(trans_),
    dt(dtheta/sigma), max_t(2.*fabs(theta_v-theta_i)/sigma), max_x(6.*sqrt(max_t)), m((int)(max_t/dt)),
    theta_array(new double[m]), flux_array(new double[m]) {} 

/* The class destructor frees the dynamically allocated memory. */
butler_volmer::~butler_volmer() {
    delete [] theta_array;
    delete [] flux_array;
}


/** Initializes adaptive spatial grid. 
* \param[in] omega expansion factor.
* \param[in] h spacing between the first two grid points. */
void butler_volmer::init_grid(double omega,double h) {
	
	// Set unequally spaced grid points 
	x.push_back(0.);
	while (x.back() < max_x) {
		x.push_back(x.back() + h);
		h *= omega;
	}
	n = x.size(); 

	// Calculate Thomas coefficients 
	alpha.assign(n, 0.);
	beta.assign(n, 0.);
	gamma.assign(n, 0.);
	for (int i=1; i<n-1; i++) {
		double dx_m = x[i] - x[i-1];
		double dx_p = x[i+1] - x[i];
		alpha[i] = -2. * dt / ( dx_m * (dx_m + dx_p) );
		gamma[i] = -2. * dt / ( dx_p * (dx_m + dx_p) );
		beta[i] = 1. - alpha[i] - gamma[i];
	}

}


/** Solves the diffusion equation, storing snapshots of the solution to a file.
 * \param[in] filename the name of the file to write to.
 * \param[in] snaps the number of snapshots to save (not including the initial snapshot). */
void butler_volmer::solve(const char* filename,int snaps) {

	// Initialize solution vector 
	C.assign(n, 1.);
	g_mod.assign(n-1, 0.);

	if (snaps > 0) {

		// Allocate memory to store solution snapshots
		double *z=new double[n*(snaps+1)];
	   	memcpy(z,C.data(),n*sizeof(double));
	   	int snapi = 1;

		// Perform timesteps 
		double theta = theta_i;
		for (int k=0; k<m; k++) {
			
			// Change potential 
			if (k < m/2) { 
				theta -= dtheta; 
			} else { 
				theta += dtheta; 
			}

			// Step the solution forward 
			step_forward(theta);

			// Compute flux 
			double flux = - (-C[2]+4.*C[1]-3.*C[0]) / (2.*(x[1]-x[0]));
			theta_array[k] = theta;
			flux_array[k] = flux;

			// Store the snapshot
			if (k == m*snapi/snaps-1) {
	        	memcpy(z+snapi*n,C.data(),n*sizeof(double));
	        	snapi += 1;
			}
		}

		// Open the output file to store the snapshots
	    FILE *fp=fopen(filename,"w");
	    if (fp==NULL) {
	        fputs("Can't open output file\n",stderr);
	        exit(1);
	    }

	    // Print the snapshots
	    for (int i=0;i<n;i++) print_line(fp,x[i],z+i,snaps);

	    // Delete snapshots and close file
	    fclose(fp);
	    delete [] z;
	
	} else {

		// Perform timesteps 
		double theta = theta_i;
		for (int k=0; k<m; k++) {
			
			// Change potential 
			if (k < m/2) { 
				theta -= dtheta; 
			} else { 
				theta += dtheta; 
			}

			// Step the solution forward 
			step_forward(theta);

			// Compute flux 
			double h = x[1]-x[0], omega = (x[2]-x[1])/h; 
			double flux = -( -C[2]+(omega+1.)*(omega+1.)*C[1]-omega*(omega+2.)*C[0] ) / ( h*omega*(omega+1.) );
			theta_array[k] = theta;
			flux_array[k] = flux;
		}
	}
}


/** Steps the solution forward in time.
 * \param[in] theta current potential. */
void butler_volmer::step_forward(double theta) {
	double fac = (x[1]-x[0]) * exp(-trans*theta) * k0;

	// Apply boundary conditions specified by Butler_Volmer kinetics 
	alpha[0] = 0.;
	beta[0] = 1. + fac * (1+exp(theta));
	gamma[0] = -1.;
	C[0] = fac * exp(theta);
	C[0] /= beta[0]; 

	// Since g_mod[0] is no longer constant, needs to update other g_mod accordingly
	g_mod[0] = gamma[0] / beta[0];
	for (int i=1; i<n-1; i++) {
		g_mod[i] = gamma[i] / (beta[i] - g_mod[i-1] * alpha[i]);
	}

	// Forward sweep 
	for (int i=1; i<n-1; i++) {
		C[i] = (C[i] - C[i-1] * alpha[i]) / (beta[i] - g_mod[i-1] * alpha[i]);
	}
	
	// Back substitution
	for (int i=n-2; i>=0; i--) {
		C[i] -= g_mod[i] * C[i+1];
	}
}


/** Prints a line of stored snapshots to a file.
 * \param[in] fp a pointer to the file to write to.
 * \param[in] xx the position in the domain corresponding to this line.
 * \param[in] zp a pointer to the first snapshot data point to print.
 * \param[in] snaps the number of snapshots (not including the starting
 *                  snapshot). */
void butler_volmer::print_line(FILE *fp,double xx,double *zp,int snaps) {
    fprintf(fp,"%g",xx);
    for(int i=0; i<=snaps; i++) fprintf(fp," %g",zp[i*n]);
    fputc('\n',fp);
}


/** Outputs cyclic voltammetry data to a file. */
void butler_volmer::output_cv(const char* filename) {
    FILE *fp=fopen(filename,"w");
    if (fp==NULL) {
        fputs("Can't open output file\n",stderr);
        exit(1);
    }
    for (int k=0; k<m; k++) fprintf(fp, "%g %g\n", theta_array[k], flux_array[k]);

    fclose(fp);
}


/** Finds postive and negative peaks in the cyclic voltammogram and prints corresponding flux and theta values, 
 *	as well as peak separation in theta. */
void butler_volmer::print_peaks() {
	double flux_pos=-100., flux_neg=100., theta_pos=0., theta_neg=0., peak_separation=0.;

	for (int i=0; i<m; i++) {
		if (flux_array[i] > flux_pos) {
			flux_pos = flux_array[i];
			theta_pos = theta_array[i];
		}
		if (flux_array[i] < flux_neg) { 
			flux_neg = flux_array[i];
			theta_neg = theta_array[i];
		}
	}

	peak_separation = theta_pos-theta_neg;
	printf("%g %g %g %g %g\n",flux_pos,flux_neg,theta_pos,theta_neg,peak_separation);
}
