#define INITCOND 2 //sets initial condition. 1 : smooth oblate gaussian, 2 : fluctuating oblate gaussian , 3 : MC-Glauber + longitudinal form , 4 : read in from file 
#define PI 3.141592654f //question bro?
#define SIGMA 1.0f //the width of the gaussian distribution in (y - eta). Generalize this to a spacetime function?
#define DIM_X 51 //number of grid points in x direction
#define DIM_Y 51 //number of grid points in y direction
#define DIM_ETA 51 //number of grid points in eta (spacetime rapidity) direction
#define DIM (DIM_X * DIM_Y * DIM_ETA) //total number of spatial grid points
#define DIM_RAP 31 //number of grid points in momentum rapidity
#define DIM_PHIP 31 //number of grid points in phi_p momentum azimuthal angle
#define DX 0.05f //spacing of grid in x direction
#define DY 0.05f //spacing of grid in y direction
#define DETA 0.05f //spacing of grid in eta direction
#define DRAP 0.1f //spacing of grid in momentum rapidity
#define DTAU 0.1f //free streaming longitudinal proper time step size
#define TAU0 0.1f //initial longitudinal proper time
#define TAU (TAU0 + DTAU) //final matching longitudinal proper time
#define EOS_TYPE 2 // 1 for conformal EOS, 2 for Wuppertal-Budhapest parameterization
