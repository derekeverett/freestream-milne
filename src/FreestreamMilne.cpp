//This file contains a wrapper class for freestream-milne

#ifndef SRC_FREESTREAMMILNE_
#define SRC_FREESTREAMMILNE_

//#include "Parameter.h"
#include "FreeStream.cpp"
#include "InitialConditions.cpp"
#include "LandauMatch.cpp"
#include "EquationOfState.cpp"
#include "HydroValidity.cpp"
#include "Memoryf.cpp"
#include "FileIO.cpp"
//#include "WriteHistory.cpp"
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#define PI 3.141592654f
#define PRINT_SCREEN 1 //turn on for program info to print to terminal
//#define HBARC 0.197326938f
#define HBARC 0.197326938 //used to convert units of input / output hydro vectors
#define TEST_INTERPOL 0 //try calculating the profile at intermediate times by interpolating between the initial profile and final (streamed) profile...

using namespace std;

class FREESTREAMMILNE {
 private:

 public:
    FREESTREAMMILNE();
    ~FREESTREAMMILNE();

    int run_freestream_milne();

    // IS THIS VARIABLE NECESSARY
    int gridSize; //the total number of grid points in x, y, and eta : used for vector memory allocation

    //support to initilialize the energy density from a vector - useful for JETSCAPE
    //note units of argument should be GeV / fm^3
    //then we convert to fm^(-4)
    void initialize_from_vector(std::vector<float>);
    std::vector<float> init_energy_density;

    //support to write final hydro variables to vectors - useful for JETSCAPE
    //note we need to convert back to GeV / fm^3 units here
    void output_to_vectors(std::vector<double>&, //e
                            std::vector<double>&, //p
                            std::vector<double>&, //ut
                            std::vector<double>&, //ux
                            std::vector<double>&, //uy
                            std::vector<double>&, //un
                            std::vector<double>&, //pitt
                            std::vector<double>&, //pitx
                            std::vector<double>&, //pity
                            std::vector<double>&, //pitn
                            std::vector<double>&, //pixx
                            std::vector<double>&, //pixy
                            std::vector<double>&, //pixn
                            std::vector<double>&, //piyy
                            std::vector<double>&, //piyn
                            std::vector<double>&, //pinn
                            std::vector<double>&); //Pi

    std::vector<double> final_energy_density;
    std::vector<double> final_pressure;
    std::vector<double> final_ut;
    std::vector<double> final_ux;
    std::vector<double> final_uy;
    std::vector<double> final_un;
    std::vector<double> final_pitt;
    std::vector<double> final_pitx;
    std::vector<double> final_pity;
    std::vector<double> final_pitn;
    std::vector<double> final_pixx;
    std::vector<double> final_pixy;
    std::vector<double> final_pixn;
    std::vector<double> final_piyy;
    std::vector<double> final_piyn;
    std::vector<double> final_pinn;
    std::vector<double> final_Pi;

};

FREESTREAMMILNE::FREESTREAMMILNE() {

}

FREESTREAMMILNE::~FREESTREAMMILNE() {
}

//use this function to initialize energy density within JETSCAPE
void FREESTREAMMILNE::initialize_from_vector(std::vector<float> energy_density_in) {
  init_energy_density = energy_density_in;
}

//use this function to return final hydro variables as vectors within JETSCAPE
void FREESTREAMMILNE::output_to_vectors(std::vector<double> &energy_density_out,
                                        std::vector<double> &pressure_out,
                                        std::vector<double> &ut_out,
                                        std::vector<double> &ux_out,
                                        std::vector<double> &uy_out,
                                        std::vector<double> &un_out,
                                        std::vector<double> &pitt_out,
                                        std::vector<double> &pitx_out,
                                        std::vector<double> &pity_out,
                                        std::vector<double> &pitn_out,
                                        std::vector<double> &pixx_out,
                                        std::vector<double> &pixy_out,
                                        std::vector<double> &pixn_out,
                                        std::vector<double> &piyy_out,
                                        std::vector<double> &piyn_out,
                                        std::vector<double> &pinn_out,
                                        std::vector<double> &Pi_out) {
  energy_density_out = final_energy_density;
  pressure_out = final_pressure;
  ut_out = final_ut;
  ux_out = final_ux;
  uy_out = final_uy;
  un_out = final_un;
  pitt_out = final_pitt;
  pitx_out = final_pitx;
  pity_out = final_pity;
  pitn_out = final_pitn;
  pixx_out = final_pixx;
  pixy_out = final_pixy;
  pixn_out = final_pixn;
  piyy_out = final_piyy;
  piyn_out = final_piyn;
  pinn_out = final_pinn;
  Pi_out = final_Pi;
}

//where the magic happens
int FREESTREAMMILNE::run_freestream_milne() {

if(PRINT_SCREEN) printf("Welcome to freestream-milne\n");

//declare parameter struct
struct parameters params;

//set default parameters in case of missing freestream_input file
params.OUTPUTFORMAT = 2;
params.BARYON = 0;
params.IC_ENERGY = 5;
params.IC_BARYON = 1;
params.ETA_WIDTH = 0.5;
params.ETA_FLAT = 0.5;
params.SIGMA = 1.0;
params.SIGMA_B = 1.0;
params.DIM_X = 101;
params.DIM_Y = 101;
params.DIM_ETA = 1;
params.DIM_RAP = 1;
params.DIM_PHIP = 51;
params.DX = 0.1;
params.DY = 0.1;
params.DETA = 0.1;
params.DRAP = 0.2;
params.DTAU = 0.5;
params.TAU0 = 0.1;
params.EOS_TYPE = 1;
params.E_FREEZE = 1.7;
params.VISCOUS_MATCHING = 1;

//read in chosen parameters from freestream_input if such a file exists
readInParameters(params);

//define some useful combinations
params.DIM = params.DIM_X * params.DIM_Y * params.DIM_ETA;
params.TAU = params.TAU0 + params.DTAU;

int DIM_X = params.DIM_X;
int DIM_Y = params.DIM_Y;
int DIM_ETA = params.DIM_ETA;

if(PRINT_SCREEN)
  {
    printf("Parameters are ...\n");
    printf("(DIM_X, DIM_Y, DIM_ETA, DIM_PHIP, DIM_RAP) = (%d, %d, %d, %d, %d)\n", params.DIM_X, params.DIM_Y, params.DIM_ETA, params.DIM_PHIP, params.DIM_RAP);
    printf("(DX, DY, DETA, DTAU) = (%.2f fm, %.2f fm, %.2f, %.2f fm/c)\n", params.DX, params.DY, params.DETA, params.DTAU);
    printf("TAU0 = %.2f fm/c\n", params.TAU0);
    printf("SIGMA = %.2f \n", params.SIGMA);
    printf("E_FREEZE = %.3f GeV / fm^3 \n", params.E_FREEZE);
    if (params.VISCOUS_MATCHING) printf("Will match to hydro including viscous part of Tmunu \n");
    else printf("Will match to hydro with ideal part of Tmunu \n");
    if (params.EOS_TYPE == 1) printf("Using EoS : Conformal \n");
    else if (params.EOS_TYPE == 2) printf("Using EoS : Wuppertal-Budhapest \n");
    else if (params.EOS_TYPE == 3) printf("Using EoS : Lattice QCD + HRG matched.\n");
  }
//allocate and initialize memory
if (PRINT_SCREEN) printf("Allocating memory\n");

//the initial energy density spatial profile
float *initialEnergyDensity = NULL;
initialEnergyDensity = (float *)calloc(params.DIM, sizeof(float));

//the initial baryon density spatial profile
float *initialChargeDensity = NULL;
if(params.BARYON) initialChargeDensity = (float *)calloc(params.DIM, sizeof(float));

//the initial density G(tilde)^(tau,tau) at time tau_0
float **density = NULL;
density = calloc2dArrayf(density, params.DIM, params.DIM_RAP); // function of x,y,eta and rapidity

//the initial density J(tilde)^(tau) at time tau_0
float **chargeDensity = NULL;
if(params.BARYON) chargeDensity = calloc2dArrayf(chargeDensity, params.DIM, params.DIM_RAP); // function of x,y,eta and rapidity

//initialize energy density

//define a lower bound on energy density for all cells to regulate numerical noise in flow velocity in dilute regions
float lower_tolerance = 1.0e-7;

if (PRINT_SCREEN) printf("setting initial conditions on energy density : ");
if (params.IC_ENERGY == 1)
{
  initializeEllipticalGauss(initialEnergyDensity, 15.0, 15.0, 15.0, params);
  if(PRINT_SCREEN) printf("Smooth Oblate Gaussian \n");
}
else if (params.IC_ENERGY == 2)
{
  initializeEllipticalMCGauss(initialEnergyDensity, 15.0, 15.0, 15.0, params);
  if(PRINT_SCREEN) printf("Fluctuating Oblate Gaussian \n");
}
else if (params.IC_ENERGY == 3)
{
  readDensityFile(initialEnergyDensity, (char *)"initial_profiles/e", params);
  if(PRINT_SCREEN) printf("Reading from energy density file in initial_profiles/ \n");
}
else if (params.IC_ENERGY == 4)
{
  readEnergyDensitySuperMCBlock(initialEnergyDensity, params);
  if(PRINT_SCREEN) printf("Reading from superMC energy density file in initial_profiles/ \n");
}
else if (params.IC_ENERGY == 5)
{
  //read in initial energy density using the initiliaze_from_vector() function
  //note that this is not safe - if one passes an empty vector it will not throw an error
  //converting units of energy density from GeV / fm^3 to fm^(-4)
  if(PRINT_SCREEN) printf("Reading energy density from initial energy density vector\n");
  //do a value copy

  //try adding a small value everywhere to regulate problems with flow velocity in dilute regions

  //TEMPORARY
  //rescale initial distribution
  float rescale = 1.0;
  //TEMPORARY
  for (int i = 0; i < params.DIM; i++) initialEnergyDensity[i] = init_energy_density[i] * rescale / (float)HBARC + lower_tolerance;
  //for (int i = 0; i < params.DIM; i++) initialEnergyDensity[i] = init_energy_density[i] / (float)HBARC;
  //just doing this here for testing - try increasing normalization of initial distribution to improve stability
  //for (int i = 0; i < params.DIM; i++) initialEnergyDensity[i] = init_energy_density[i];
}
else
{
  printf("Not a valid initial Condition... Goodbye\n");
  return 0;
}

if (params.BARYON)
{
  //initialize baryon density
  if (PRINT_SCREEN) printf("setting initial conditions on baryon density : ");
  if (params.IC_BARYON == 1)
  {
    initializeEllipticalGauss(initialChargeDensity, 2.0, 3.0, 1.0, params);
    if (PRINT_SCREEN) printf("Smooth Oblate Gaussian \n");
  }
  else if (params.IC_BARYON == 2)
  {
    initializeEllipticalMCGauss(initialChargeDensity, 2.0, 3.0, 1.0, params);
    if (PRINT_SCREEN) printf("Fluctuating Oblate Gaussian \n");
  }
  else if (params.IC_BARYON == 3)
  {
    readDensityFile(initialChargeDensity, (char *)"initial_profiles/nB", params);
    if (PRINT_SCREEN) printf("Reading from baryon density file in initial_profiles/ \n");
  }
  else
  {
    printf("Not a valid initial Condition... Goodbye\n");
    return 0;
  }
}

//write initial energy density and baryon density to file
writeScalarToFile(initialEnergyDensity, (char *)"initial_e", params);
if (params.BARYON) writeScalarToFile(initialChargeDensity, (char *)"initial_nB", params);
writeScalarToFileProjection(initialEnergyDensity, (char *)"initial_e_projection", params);
if (params.BARYON) writeScalarToFileProjection(initialChargeDensity, (char *)"initial_nB_projection", params);

/////////////////////////////BEGIN TESTING FOR JETSCAPE//////////////////////////////
//make a toy plot of 1/tau * initial energy density to compare 2+1D freestreaming with only longitudinal (bjorken) dilution
float *scaledEnergyDensity = NULL;
if (TEST_INTERPOL)
{
  printf("Calculating 1 / tau scaled profile for testing \n");
  scaledEnergyDensity = (float *)calloc(params.DIM, sizeof(float));
  for (int is = 0; is < params.DIM; is++) scaledEnergyDensity[is] = initialEnergyDensity[is] * (params.TAU0 / params.TAU);
  writeScalarToFile(scaledEnergyDensity, (char *)"scaled_e", params);
  writeScalarToFileProjection(scaledEnergyDensity, (char *)"scaled_e_projection", params);
}
/////////////////////////////END TESTING FOR JETSCAPE//////////////////////////////


//calculate total energy to check convergence
float totalEnergy = 0.0;
for (int is = 0; is < params.DIM; is++) totalEnergy += initialEnergyDensity[is];
if (params.DIM_ETA > 1) totalEnergy *= (params.TAU0 * params.DX * params.DY * params.DETA);
else totalEnergy *= (params.DX * params.DY);
printf("Total energy before streaming : %f \n", totalEnergy);

//convert the energy density profile into the initial density profile to be streamed and free memory
convertInitialDensity(initialEnergyDensity, density, params);
if (!TEST_INTERPOL) free(initialEnergyDensity);
//convert the baryon density profile into the initial baryon density profile to be streamed and free memory
if (params.BARYON) convertInitialChargeDensity(initialChargeDensity, chargeDensity, params);
if (params.BARYON) free(initialChargeDensity);

//the shifted energy density profile G^(tau,tau) at time tau
float ***shiftedDensity = NULL;
shiftedDensity = calloc3dArrayf(shiftedDensity, params.DIM, params.DIM_RAP, params.DIM_PHIP);

//the shifted baryon density profile J^(tau) at time tau
float ***shiftedChargeDensity = NULL;
if(params.BARYON) shiftedChargeDensity = calloc3dArrayf(shiftedChargeDensity, params.DIM, params.DIM_RAP, params.DIM_PHIP);

//perform the free streaming time-update step and free up memory
//pretabulate trig and hypertrig functions before this step to save time?
if (PRINT_SCREEN) printf("performing the free streaming\n");

double sec = 0.0;
#ifdef _OPENMP
sec = omp_get_wtime();
#endif

///////////  BEGIN LOOP OVER TIME STEPS HERE ////////////////////////
////////// MOVE DECLARATIONS AND ALLOCATION OUTSIDE LOOP? //////////
freeStream(density, shiftedDensity, params);

//#pragma acc update host(shiftedDensity)
free2dArrayf(density, params.DIM);
if (params.BARYON) freeStream(chargeDensity, shiftedChargeDensity, params);
//#if(params.BARYON) pragma acc update host(shiftedChargeDensity)
if (params.BARYON) free2dArrayf(chargeDensity, params.DIM);

#ifdef _OPENMP
sec = omp_get_wtime() - sec;
#endif
if (PRINT_SCREEN) printf("Free streaming took %f seconds\n", sec);

//Landau matching to find the components of energy-momentum tensor
if (PRINT_SCREEN) printf("Landau matching to find hydrodynamic variables\n");

//the ten independent components of the stress tensor
float **stressTensor = NULL;
stressTensor = calloc2dArrayf(stressTensor, 10, params.DIM);

//the four independent components of baryon current four-vector
float **baryonCurrent = NULL;
if(params.BARYON) baryonCurrent = calloc2dArrayf(baryonCurrent, 4, params.DIM);

//a table containing 10 rows for 10 independent combinations of p_(mu)p_(nu)
//hypertrig table depends on TAU, so need to keep this inside loop
float ****hypertrigTable = NULL;
hypertrigTable = calloc4dArrayf(hypertrigTable, 10, params.DIM_RAP, params.DIM_PHIP, params.DIM_ETA); //depends on eta because we have function of eta - y

if (PRINT_SCREEN) printf("calculating hypertrig table\n");
calculateHypertrigTable(hypertrigTable, params);

//calculate the ten independent components of the stress tensor by integrating over rapidity and phi_p
if (PRINT_SCREEN) printf("calculating independent components of stress tensor\n");
#ifdef _OPENMP
sec = omp_get_wtime();
#endif
calculateStressTensor(stressTensor, shiftedDensity, hypertrigTable, params);
free3dArrayf(shiftedDensity, params.DIM, params.DIM_RAP);
#ifdef _OPENMP
sec = omp_get_wtime() - sec;
#endif
if (PRINT_SCREEN) printf("calculating stress tensor took %f seconds\n", sec);

if (params.BARYON)
{
  //calculate the four independent components of the baryon current by integrating over rapidity and phi_p
  if (PRINT_SCREEN) printf("calculating independent components of baryon current\n");
  #ifdef _OPENMP
  sec = omp_get_wtime();
  #endif
  calculateBaryonCurrent(baryonCurrent, shiftedChargeDensity, hypertrigTable, params);
  free3dArrayf(shiftedChargeDensity, params.DIM, params.DIM_RAP);
  #ifdef _OPENMP
  sec = omp_get_wtime() - sec;
  #endif
  if (PRINT_SCREEN) printf("calculating baryon current took %f seconds\n", sec);
}

//done with hypertrig table as well
free4dArrayf(hypertrigTable, 10, params.DIM_RAP, params.DIM_PHIP);

//variables to store the hydrodynamic variables after the Landau matching is performed
//the energy density
float *energyDensity = NULL;
energyDensity = (float *)calloc(params.DIM, sizeof(float));

//the baryon density
float *baryonDensity = NULL;
if(params.BARYON) baryonDensity = (float *)calloc(params.DIM, sizeof(float));

//the flow velocity
float **flowVelocity = NULL;
flowVelocity = calloc2dArrayf(flowVelocity, 4, params.DIM);

//the pressure
float *pressure = NULL;
pressure = (float *)calloc(params.DIM, sizeof(float));

//the bulk pressure Pi
float *bulkPressure = NULL;
bulkPressure = (float *)calloc(params.DIM, sizeof(float));

//the shear stress tensor
float **shearTensor = NULL;
shearTensor = calloc2dArrayf(shearTensor, 10, params.DIM); //calculate 10 components, can check tracelessness/orthogonality for accuracy

//the baryon diffusion current vector
float **baryonDiffusion = NULL;
if(params.BARYON) baryonDiffusion = calloc2dArrayf(baryonDiffusion, 4, params.DIM);

//solve the eigenvalue problem for the energy density and flow velocity
if (PRINT_SCREEN) printf("solving eigenvalue problem for energy density and flow velocity\n");
#ifdef _OPENMP
sec = omp_get_wtime();
#endif
solveEigenSystem(stressTensor, energyDensity, flowVelocity, params);
#ifdef _OPENMP
sec = omp_get_wtime() - sec;
#endif
if (PRINT_SCREEN) printf("solving eigenvalue problem took %f seconds\n", sec);

if (params.BARYON)
{
  //calculate baryon density and diffusion current
  if (PRINT_SCREEN) printf("calculating baryon density and diffusion current \n");
  #ifdef _OPENMP
  sec = omp_get_wtime();
  #endif
  calculateBaryonDensity(baryonDensity, baryonCurrent, flowVelocity, params);
  calculateBaryonDiffusion(baryonDiffusion, baryonCurrent, baryonDensity, flowVelocity, params);
  #ifdef _OPENMP
  sec = omp_get_wtime() - sec;
  #endif
  if (PRINT_SCREEN) printf("calculating baryon density and diffusion current took %f seconds\n", sec);
}

calculatePressure(energyDensity, baryonDensity, pressure, params);

//viscous currents are nonzero if VISCOUS_MATCHING is on
if (params.VISCOUS_MATCHING)
{
  calculateBulkPressure(stressTensor, energyDensity, pressure, bulkPressure, params);
  calculateShearViscTensor(stressTensor, energyDensity, flowVelocity, pressure, bulkPressure, shearTensor, params);
}


/////////////////////////////BEGIN TESTING FOR JETSCAPE//////////////////////////////
if (TEST_INTERPOL)
{
  printf("approximating energy density profile at intermed. times by interpolating between initial and final profiles \n");
  float TAU = params.TAU;
  float TAU0 = params.TAU0;
  float DTAU = params.DTAU;
  float tau_i = TAU0 + (DTAU / 2.0); //some intermediate time
  float c_1 = (TAU0 / tau_i);
  float c_2 = (tau_i - TAU0) / DTAU / tau_i;
  for (int is = 0; is < params.DIM; is++) scaledEnergyDensity[is] = c_1 * initialEnergyDensity[is] + c_2 * ((TAU * energyDensity[is] - TAU0 * initialEnergyDensity[is]));
  writeScalarToFile(scaledEnergyDensity, (char *)"tau_interpolated_e", params);
  writeScalarToFileProjection(scaledEnergyDensity, (char *)"tau_interpolated_e_projection", params);
  free(scaledEnergyDensity);
}
/////////////////////////////END TESTING FOR JETSCAPE//////////////////////////////

//write the next entry in hdf5 spacetime evolution file
//WriteHistoryJETSCAPE();

///////////  END LOOP OVER TIME STEPS HERE ////////////////////////

float totalEnergyAfter = 0.0;
for (int is = 0; is < params.DIM; is++) totalEnergyAfter += stressTensor[0][is];
if (params.DIM_ETA > 1) totalEnergyAfter *= (params.TAU * params.DX * params.DY * params.DETA);
else totalEnergyAfter *= (params.TAU * params.DX * params.DY);
printf("Total energy after streaming : %f \n", totalEnergyAfter);

//check which fraction of total energy lies within freezeout surface, which lies in 'corona'
float totalEnergyInsideHypersurf = 0.0;
for (int is = 0; is < params.DIM; is++)
{
  if ( (energyDensity[is] * HBARC) > params.E_FREEZE) totalEnergyInsideHypersurf += energyDensity[is];
}
if (params.DIM_ETA > 1) totalEnergyInsideHypersurf *= (params.TAU * params.DX * params.DY * params.DETA);
else totalEnergyInsideHypersurf *= (params.TAU * params.DX * params.DY);
printf("Fraction of energy contained in Freezeout Hypersurface : %f \n", totalEnergyInsideHypersurf / totalEnergyAfter);

//////////////////////////////////HYDRO VALIDITY//////////////////////////////////
//bulk inv reynolds #
float *R_Pi_Inv = NULL;
R_Pi_Inv = (float *)calloc(params.DIM, sizeof(float));
//shear inv reynolds #
float *R_pimunu_Inv = NULL;
R_pimunu_Inv = (float *)calloc(params.DIM, sizeof(float));
calculateBulkInvReynolds(pressure, bulkPressure, R_Pi_Inv, params);
calculateShearInvReynolds(energyDensity, pressure, shearTensor, R_pimunu_Inv, params);
writeScalarToFileProjection(R_Pi_Inv, (char *)"R_Pi_Inv_projection", params);
writeScalarToFileProjection(R_pimunu_Inv, (char *)"R_pimunu_Inv_projection", params);
int ctr = (DIM_Y * DIM_ETA * ((DIM_X - 1) / 2)) + (DIM_ETA * ((DIM_Y - 1) / 2)) + ((DIM_ETA - 1) / 2);
printf("R_Pi_Inv at center : %f \n", R_Pi_Inv[ctr]);
printf("R_pimunu_Inv at center : %f \n", R_pimunu_Inv[ctr]);
free(R_Pi_Inv);
free(R_pimunu_Inv);


//check transversality and tracelesness
//components of pi^munu u_mu
float *pi_dot_u_tau = NULL;
pi_dot_u_tau = (float *)calloc(params.DIM, sizeof(float));
float *pi_dot_u_x = NULL;
pi_dot_u_x = (float *)calloc(params.DIM, sizeof(float));
float *pi_dot_u_y = NULL;
pi_dot_u_y = (float *)calloc(params.DIM, sizeof(float));
float *pi_dot_u_eta = NULL;
pi_dot_u_eta = (float *)calloc(params.DIM, sizeof(float));

calculate_pi_dot_u(flowVelocity, shearTensor, pi_dot_u_tau, pi_dot_u_x, pi_dot_u_y, pi_dot_u_eta, params);
writeScalarToFileProjection(pi_dot_u_tau, (char *)"pi_dot_u_tau_projection", params);
writeScalarToFileProjection(pi_dot_u_x, (char *)"pi_dot_u_x_projection", params);
writeScalarToFileProjection(pi_dot_u_y, (char *)"pi_dot_u_y_projection", params);
writeScalarToFileProjection(pi_dot_u_eta, (char *)"pi_dot_u_eta_projection", params);

free(pi_dot_u_tau);
free(pi_dot_u_x);
free(pi_dot_u_y);
free(pi_dot_u_eta);

//trace of pi^munu , pi^mu_mu
float *pi_mu_mu = NULL;
pi_mu_mu = (float *)calloc(params.DIM, sizeof(float));

calculate_pi_mu_mu(shearTensor, pi_mu_mu, params);
writeScalarToFileProjection(pi_mu_mu, (char *)"pi_mu_mu_projection", params);
free(pi_mu_mu);

//////////////////////////////////HYDRO VALIDITY//////////////////////////////////


if (PRINT_SCREEN) printf("writing hydro variables\n");

writeScalarToFile(energyDensity, (char *)"e", params);
writeScalarToFile(pressure, (char *)"p", params);
writeScalarToFile(bulkPressure, (char *)"bulk_PI", params);
writeScalarToFileProjection(energyDensity, (char *)"e_projection", params);
writeScalarToFileProjection(pressure, (char *)"p_projection", params);
writeScalarToFileProjection(bulkPressure, (char *)"bulk_PI_projection", params);

writeVectorToFile(flowVelocity, (char *)"u_tau", 0, params);
writeVectorToFile(flowVelocity, (char *)"u_x", 1, params);
writeVectorToFile(flowVelocity, (char *)"u_y", 2,params);
writeVectorToFile(flowVelocity, (char *)"u_eta", 3,params);

writeVectorToFileProjection(flowVelocity, (char *)"u_tau_projection", 0,params);
writeVectorToFileProjection(flowVelocity, (char *)"u_x_projection", 1,params);
writeVectorToFileProjection(flowVelocity, (char *)"u_y_projection", 2,params);
writeVectorToFileProjection(flowVelocity, (char *)"u_eta_projection", 3,params);


writeVectorToFile(shearTensor, (char *)"pi_tau_tau", 0,params);
writeVectorToFile(shearTensor, (char *)"pi_tau_x", 1,params);
writeVectorToFile(shearTensor, (char *)"pi_tau_y", 2,params);
writeVectorToFile(shearTensor, (char *)"pi_tau_eta", 3,params);
writeVectorToFile(shearTensor, (char *)"pi_x_x", 4,params);
writeVectorToFile(shearTensor, (char *)"pi_x_y", 5,params);
writeVectorToFile(shearTensor, (char *)"pi_x_eta", 6,params);
writeVectorToFile(shearTensor, (char *)"pi_y_y", 7,params);
writeVectorToFile(shearTensor, (char *)"pi_y_eta", 8,params);
writeVectorToFile(shearTensor, (char *)"pi_eta_eta", 9,params);

writeVectorToFileProjection(shearTensor, (char *)"pi_tau_tau_projection", 0,params);
writeVectorToFileProjection(shearTensor, (char *)"pi_tau_x_projection", 1,params);
writeVectorToFileProjection(shearTensor, (char *)"pi_tau_y_projection", 2,params);
writeVectorToFileProjection(shearTensor, (char *)"pi_tau_eta_projection", 3,params);
writeVectorToFileProjection(shearTensor, (char *)"pi_x_x_projection", 4,params);
writeVectorToFileProjection(shearTensor, (char *)"pi_x_y_projection", 5,params);
writeVectorToFileProjection(shearTensor, (char *)"pi_x_eta_projection", 6,params);
writeVectorToFileProjection(shearTensor, (char *)"pi_y_y_projection", 7,params);
writeVectorToFileProjection(shearTensor, (char *)"pi_y_eta_projection", 8,params);
writeVectorToFileProjection(shearTensor, (char *)"pi_eta_eta_projection", 9,params);

if (params.BARYON)
{
  writeScalarToFile(baryonDensity, (char *)"nB",params);
  writeScalarToFileProjection(baryonDensity, (char *)"nB_projection",params);
  writeVectorToFile(baryonDiffusion, (char *)"V_x", 1,params);
  writeVectorToFile(baryonDiffusion, (char *)"V_y", 2,params);
  writeVectorToFile(baryonDiffusion, (char *)"V_eta", 3,params);
}


//support for JETSCAPE - write hydro variables to vectors
final_energy_density.resize(params.DIM);
final_pressure.resize(params.DIM);
final_ut.resize(params.DIM);
final_ux.resize(params.DIM);
final_uy.resize(params.DIM);
final_un.resize(params.DIM);
final_pitt.resize(params.DIM);
final_pitx.resize(params.DIM);
final_pity.resize(params.DIM);
final_pitn.resize(params.DIM);
final_pixx.resize(params.DIM);
final_pixy.resize(params.DIM);
final_pixn.resize(params.DIM);
final_piyy.resize(params.DIM);
final_piyn.resize(params.DIM);
final_pinn.resize(params.DIM);
final_Pi.resize(params.DIM);

if ( (params.OUTPUTFORMAT == 2) || (params.OUTPUTFORMAT == 3) )
{
  for (int is = 0; is < params.DIM; is++)
  {
    //converting back to GeV / fm^3 for use in JETSCAPE
    final_energy_density[is] = (double)energyDensity[is] * HBARC;
    final_pressure[is] = (double)pressure[is] * HBARC;
    final_ut[is] = (double)flowVelocity[0][is];
    final_ux[is] = (double)flowVelocity[1][is];
    final_uy[is] = (double)flowVelocity[2][is];
    final_un[is] = (double)flowVelocity[3][is];
    final_pitt[is] = (double)shearTensor[0][is] * HBARC;
    final_pitx[is] = (double)shearTensor[1][is] * HBARC;
    final_pity[is] = (double)shearTensor[2][is] * HBARC;
    final_pitn[is] = (double)shearTensor[3][is] * HBARC;
    final_pixx[is] = (double)shearTensor[4][is] * HBARC;
    final_pixy[is] = (double)shearTensor[5][is] * HBARC;
    final_pixn[is] = (double)shearTensor[6][is] * HBARC;
    final_piyy[is] = (double)shearTensor[7][is] * HBARC;
    final_piyn[is] = (double)shearTensor[8][is] * HBARC;
    final_pinn[is] = (double)shearTensor[9][is] * HBARC;
    final_Pi[is] = (double)bulkPressure[is] * HBARC;
  }
}

//free the memory
free2dArrayf(stressTensor, 10);
free(energyDensity);
free2dArrayf(flowVelocity, 4);
free(pressure);
free(bulkPressure);
free2dArrayf(shearTensor, 10);

if (params.BARYON)
{
  free2dArrayf(baryonCurrent, 4);
  free(baryonDensity);
  free2dArrayf(baryonDiffusion, 4);
}

if (PRINT_SCREEN) printf("Done... Goodbye!\n");

//change this to return a different int status if something goes wrong?
return 0;
}

#endif  // SRC_FREESTREAMMILNE_
