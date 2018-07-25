#ifndef BUTLER_VOLMER_HH
#define BUTLER_VOLMER_HH

#include <cstdio>
#include <vector> // Header for vectors


/** The solver for the butler_volmer model of the electrochemical CV process. 
 *  Backward Euler method with adaptive spatial grid and Thomas algorithm are used. 
 */
class butler_volmer {
    public:
        /** Initial potential. */
        const double theta_i;
        /** Vertex potential. */
        const double theta_v;
        /** Scan rate. */
        const double sigma;
        /** Size of potential step. */
        const double dtheta;
        /** Standard heterogeneous rate. */
        const double k0;
        /** Charge transfer coefficient. */
        const double trans;
        /** Size of timestep. */
        const double dt;
        /** Time duration to simulate. */
        const double max_t;
        /** Upper spatial boundary (lower spatial boundary is zero by default). */
        const double max_x;
        /** Number of timesteps. */
        const int m;
        /** Number of spacesteps. */
        int n;
        /** Spatial grid positions. */
        std::vector<double> x;
        /** Left, center and right coefficients for finite difference second-order derivative. */
        std::vector<double> alpha;
        std::vector<double> beta;
        std::vector<double> gamma;
        /** Thomas algorithm coefficients. */
        std::vector<double> g_mod;
        /** Solution vector. */
        std::vector<double> C;
        /** History of potential. */
        double* theta_array;
        /** History of flux. */
        double* flux_array;
        butler_volmer(double theta_i_,double theta_v_,double sigma_,double dtheta_,double k0_,double trans_);
        ~butler_volmer();
        void init_grid(double omega,double h);
        void solve(const char* filename,int snaps);
        void step_forward(double theta);
        void print_line(FILE *fp,double xx,double *zp,int snaps);
        void output_cv(const char* filename);
        void print_peaks(); 
};

#endif