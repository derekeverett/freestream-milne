// freestream-milne

//#include "Parameter.h"
#include "FreeStream.cpp"
#include "InitialConditions.cpp"
#include "LandauMatch.cpp"
#include "EquationOfState.cpp"
#include "Memory.cpp"
#include "FileIO.cpp"
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#define PI 3.141592654f
#define PRINT_SCREEN 1 //turn on for program info to print to terminal
#define TEST_INTERPOL 0 //try calculating the profile at intermediate times by interpolating between the initial profile and final (streamed) profile...

int main(void)
{
  if(PRINT_SCREEN) printf("Welcome to freestream-milne\n");

  //declare parameter struct
  struct parameters params;

  //set default parameters in case of missing freestream_input file
  params.OUTPUTFORMAT = 1;
  params.BARYON = 0;
  params.IC_ENERGY = 1;
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

  //read in chosen parameters from freestream_input if such a file exists
  readInParameters(params);

  //define some useful combinations
  params.DIM = params.DIM_X * params.DIM_Y * params.DIM_ETA;
  params.TAU = params.TAU0 + params.DTAU;

  if(PRINT_SCREEN)
    {
      printf("Parameters are ...\n");
      printf("(DIM_X, DIM_Y, DIM_ETA) = (%d, %d, %d)\n", params.DIM_X, params.DIM_Y, params.DIM_ETA);
      printf("(DX, DY, DETA, DTAU) = (%.2f, %.2f, %.2f, %.2f)\n", params.DX, params.DY, params.DETA, params.DTAU);
      if (REGULATE) printf("Regulating flow velocity in dilute region for smoothness\n");
      if (params.EOS_TYPE == 1) printf("Using EoS : Conformal \n");
      else if (params.EOS_TYPE == 2) printf("Using EoS : Wuppertal-Budhapest \n");
      else if (params.EOS_TYPE == 3) printf("Using EoS : Lattice QCD + HRG matched.\n");

    }
  //allocate and initialize memory
  if (PRINT_SCREEN) printf("Allocating memory\n");
  /*
  float ***eqnOfStateTable;
  if(params.EOS_TYPE == 3)
  { //table is regularly spaced in mu_B and T, with 450 values each
    eqnOfStateTable = calloc2dArray(eqnOfStateTable, 450, 450);
  }
  */
  //the initial energy density spatial profile


  float *initialEnergyDensity = NULL;
  initialEnergyDensity = (float *)calloc(params.DIM, sizeof(float));
  //initialEnergyDensity = (float *)malloc(params.DIM * sizeof(float));

  //the initial baryon density spatial profile
  float *initialChargeDensity = NULL;
  if(params.BARYON) initialChargeDensity = (float *)calloc(params.DIM, sizeof(float));

  //the initial density G(tilde)^(tau,tau) at time tau_0
  float **density = NULL;
  density = calloc2dArray(density, params.DIM, params.DIM_RAP); // function of x,y,eta and rapidity

  //the initial density J(tilde)^(tau) at time tau_0
  float **chargeDensity = NULL;
  if(params.BARYON) chargeDensity = calloc2dArray(chargeDensity, params.DIM, params.DIM_RAP); // function of x,y,eta and rapidity

  //initialize energy density
  if (PRINT_SCREEN) printf("setting initial conditions on energy density : ");
  if (params.IC_ENERGY == 1)
  {
    initializeEllipticalGauss(initialEnergyDensity, 10.0, 10.0, 3.0, params);
    if(PRINT_SCREEN) printf("Smooth Oblate Gaussian \n");
  }
  else if (params.IC_ENERGY == 2)
  {
    initializeEllipticalMCGauss(initialEnergyDensity, 10.0, 10.0, 3.0, params);
    if(PRINT_SCREEN) printf("Fluctuating Oblate Gaussian \n");
  }
  else if (params.IC_ENERGY == 3)
  {
    readDensityFile(initialEnergyDensity, "initial_profiles/e", params);
    if(PRINT_SCREEN) printf("Reading from energy density file in initial_profiles/ \n");
  }
  else if (params.IC_ENERGY == 4)
  {
    readEnergyDensitySuperMCBlock(initialEnergyDensity, params);
    if(PRINT_SCREEN) printf("Reading from superMC energy density file in initial_profiles/ \n");
  }
  else if (params.IC_ENERGY == 5)
  {
    //read in initial energy density within JETSCAPE framework
    if(PRINT_SCREEN) printf("Reading energy density from JETSCAPE pointer\n");
  }
  else if (params.IC_ENERGY == 6)
  {
    readEnergyDensityTRENTOBlock(initialEnergyDensity, params);
    if(PRINT_SCREEN) printf("Reading energy density from 2D TRENTO block profile in initial_profiles/\n");
  }
  else if (params.IC_ENERGY == 7)
  {
    initialize2Gaussians(initialEnergyDensity, 3.0, 3.0, 5.0, params);
    if(PRINT_SCREEN) printf("Two Oblate Gaussians \n");
  }
  else if (params.IC_ENERGY == 8)
  {
    readEnergyDensityTRENTO3DBlock(initialEnergyDensity, params);
    if(PRINT_SCREEN) printf("Reading energy density from 3D TRENTO block profile in initial_profiles/\n");
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
      readDensityFile(initialChargeDensity, "initial_profiles/nB", params);
      if (PRINT_SCREEN) printf("Reading from baryon density file in initial_profiles/ \n");
    }
    else
    {
      printf("Not a valid initial Condition... Goodbye\n");
      return 0;
    }
  }

  //write initial energy density and baryon density to file
  writeScalarToFile(initialEnergyDensity, "initial_e", params);
  if (params.BARYON) writeScalarToFile(initialChargeDensity, "initial_nB", params);
  writeScalarToFileProjection(initialEnergyDensity, "initial_e_projection", params);
  if (params.BARYON) writeScalarToFileProjection(initialChargeDensity, "initial_nB_projection", params);

  /////////////////////////////BEGIN TESTING FOR JETSCAPE//////////////////////////////
  //make a toy plot of 1/tau * initial energy density to compare 2+1D freestreaming with only longitudinal (bjorken) dilution
  float *scaledEnergyDensity = NULL;
  if (TEST_INTERPOL)
  {
    printf("Calculating 1 / tau scaled profile for testing \n");
    scaledEnergyDensity = (float *)calloc(params.DIM, sizeof(float));
    for (int is = 0; is < params.DIM; is++) scaledEnergyDensity[is] = initialEnergyDensity[is] * (params.TAU0 / params.TAU);
    writeScalarToFile(scaledEnergyDensity, "scaled_e", params);
    writeScalarToFileProjection(scaledEnergyDensity, "scaled_e_projection", params);
  }
  /////////////////////////////END TESTING FOR JETSCAPE//////////////////////////////

  //test regulating the initial profile in dilute regions
  /*
  for (int is = 0; is < params.DIM; is++)
    {
      if ( initialEnergyDensity[is] < 1.0e-10 ) initialEnergyDensity[is] = 0.0;
    }
  */
  //test regulating initial profile in dilute regions


  //convert the energy density profile into the initial density profile to be streamed and free memory
  convertInitialDensity(initialEnergyDensity, density, params);
  if (!TEST_INTERPOL) free(initialEnergyDensity);
  //convert the baryon density profile into the initial baryon density profile to be streamed and free memory
  if (params.BARYON) convertInitialChargeDensity(initialChargeDensity, chargeDensity, params);
  if (params.BARYON) free(initialChargeDensity);

  //the shifted energy density profile G^(tau,tau) at time tau
  float ***shiftedDensity = NULL;
  shiftedDensity = calloc3dArray(shiftedDensity, params.DIM, params.DIM_RAP, params.DIM_PHIP);

  //the shifted baryon density profile J^(tau) at time tau
  float ***shiftedChargeDensity = NULL;
  if(params.BARYON) shiftedChargeDensity = calloc3dArray(shiftedChargeDensity, params.DIM, params.DIM_RAP, params.DIM_PHIP);

  //perform the free streaming time-update step and free up memory
  //pretabulate trig and hypertrig functions before this step to save time?
  if (PRINT_SCREEN) printf("performing the free streaming\n");
  //copy initial and shifted density arrays to GPU
  //#pragma acc data copy(density[:params.DIM][:params.DIM_RAP]), copy(shiftedDensity[:params.DIM][:params.DIM_RAP][:params.DIM_PHIP]) //copy energy density arrays
  //#if(params.BARYON) pragma acc data copy(chargeDensity), copy(shiftedChargeDensity)  //copy baryon density arrays

  double sec = 0.0;
  #ifdef _OPENMP
  sec = omp_get_wtime();
  #endif
  freeStream(density, shiftedDensity, params);

  //#pragma acc update host(shiftedDensity)
  free2dArray(density, params.DIM);
  if (params.BARYON) freeStream(chargeDensity, shiftedChargeDensity, params);
  //#if(params.BARYON) pragma acc update host(shiftedChargeDensity)
  if (params.BARYON) free2dArray(chargeDensity, params.DIM);

  #ifdef _OPENMP
  sec = omp_get_wtime() - sec;
  #endif
  if (PRINT_SCREEN) printf("Free streaming took %f seconds\n", sec);

  //Landau matching to find the components of energy-momentum tensor
  if (PRINT_SCREEN) printf("Landau matching to find hydrodynamic variables\n");

  //the ten independent components of the stress tensor
  float **stressTensor = NULL;
  stressTensor = calloc2dArray(stressTensor, 10, params.DIM);

  //the four independent components of baryon current four-vector
  float **baryonCurrent = NULL;
  if(params.BARYON) baryonCurrent = calloc2dArray(baryonCurrent, 4, params.DIM);

  //a table containing 10 rows for 10 independent combinations of p_(mu)p_(nu)
  float ****hypertrigTable = NULL;
  hypertrigTable = calloc4dArray(hypertrigTable, 10, params.DIM_RAP, params.DIM_PHIP, params.DIM_ETA); //depends on eta because we have function of eta - y

  if (PRINT_SCREEN) printf("calculating hypertrig table\n");
  #ifdef _OPENMP
  sec = omp_get_wtime();
  #endif
  calculateHypertrigTable(hypertrigTable, params);
  #ifdef _OPENMP
  sec = omp_get_wtime() - sec;
  #endif
  if (PRINT_SCREEN) printf("calculating trig table took %f seconds\n", sec);

  //calculate the ten independent components of the stress tensor by integrating over rapidity and phi_p
  if (PRINT_SCREEN) printf("calculating independent components of stress tensor\n");
  #ifdef _OPENMP
  sec = omp_get_wtime();
  #endif
  calculateStressTensor(stressTensor, shiftedDensity, hypertrigTable, params);
  free3dArray(shiftedDensity, params.DIM, params.DIM_RAP);
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
    free3dArray(shiftedChargeDensity, params.DIM, params.DIM_RAP);
    #ifdef _OPENMP
    sec = omp_get_wtime() - sec;
    #endif
    if (PRINT_SCREEN) printf("calculating baryon current took %f seconds\n", sec);
  }

  //done with hypertrig table as well
  free4dArray(hypertrigTable, 10, params.DIM_RAP, params.DIM_PHIP);

  //variables to store the hydrodynamic variables after the Landau matching is performed
  //the energy density
  float *energyDensity = NULL;
  energyDensity = (float *)calloc(params.DIM, sizeof(float));

  //the baryon density
  float *baryonDensity = NULL;
  if(params.BARYON) baryonDensity = (float *)calloc(params.DIM, sizeof(float));

  //the flow velocity
  float **flowVelocity = NULL;
  flowVelocity = calloc2dArray(flowVelocity, 4, params.DIM);

  //the pressure
  float *pressure = NULL;
  pressure = (float *)calloc(params.DIM, sizeof(float));

  //the bulk pressure Pi
  float *bulkPressure = NULL;
  bulkPressure = (float *)calloc(params.DIM, sizeof(float));

  //the shear stress tensor
  float **shearTensor = NULL;
  shearTensor = calloc2dArray(shearTensor, 10, params.DIM); //calculate 10 components, can check tracelessness/orthogonality for accuracy

  //the baryon diffusion current vector
  float **baryonDiffusion = NULL;
  if(params.BARYON) baryonDiffusion = calloc2dArray(baryonDiffusion, 4, params.DIM);

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
  calculateBulkPressure(stressTensor, energyDensity, pressure, bulkPressure, params);
  calculateShearViscTensor(stressTensor, energyDensity, flowVelocity, pressure, bulkPressure, shearTensor, params);

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
    writeScalarToFile(scaledEnergyDensity, "tau_interpolated_e", params);
    writeScalarToFileProjection(scaledEnergyDensity, "tau_interpolated_e_projection", params);
  }
  /////////////////////////////END TESTING FOR JETSCAPE//////////////////////////////

  if (PRINT_SCREEN) printf("writing hydro variables to file\n");


  writeScalarToFile(energyDensity, "e", params);
  writeScalarToFile(pressure, "p", params);
  writeScalarToFile(bulkPressure, "bulk", params);
  writeScalarToFileProjection(energyDensity, "e_projection", params);
  writeScalarToFileProjection(pressure, "p_projection", params);
  writeScalarToFileProjection(bulkPressure, "bulk_projection", params);

  writeVectorToFile(flowVelocity, "ut", 0, params);
  writeVectorToFile(flowVelocity, "ux", 1, params);
  writeVectorToFile(flowVelocity, "uy", 2, params);
  writeVectorToFile(flowVelocity, "un", 3, params);

  writeVectorToFileProjection(flowVelocity, "ut_projection", 0, params);
  writeVectorToFileProjection(flowVelocity, "ux_projection", 1, params);
  writeVectorToFileProjection(flowVelocity, "uy_projection", 2, params);
  writeVectorToFileProjection(flowVelocity, "un_projection", 3, params);


  writeVectorToFile(shearTensor, "pitt", 0, params);
  writeVectorToFile(shearTensor, "pitx", 1, params);
  writeVectorToFile(shearTensor, "pity", 2, params);
  writeVectorToFile(shearTensor, "pitn", 3, params);
  writeVectorToFile(shearTensor, "pixx", 4, params);
  writeVectorToFile(shearTensor, "pixy", 5, params);
  writeVectorToFile(shearTensor, "pixn", 6, params);
  writeVectorToFile(shearTensor, "piyy", 7, params);
  writeVectorToFile(shearTensor, "piyn", 8, params);
  writeVectorToFile(shearTensor, "pinn", 9, params);

  writeVectorToFileProjection(shearTensor, "pitt_projection", 0, params);
  writeVectorToFileProjection(shearTensor, "pitx_projection", 1, params);
  writeVectorToFileProjection(shearTensor, "pity_projection", 2, params);
  writeVectorToFileProjection(shearTensor, "pitn_projection", 3, params);
  writeVectorToFileProjection(shearTensor, "pixx_projection", 4, params);
  writeVectorToFileProjection(shearTensor, "pixy_projection", 5, params);
  writeVectorToFileProjection(shearTensor, "pixn_projection", 6, params);
  writeVectorToFileProjection(shearTensor, "piyy_projection", 7, params);
  writeVectorToFileProjection(shearTensor, "piyn_projection", 8, params);
  writeVectorToFileProjection(shearTensor, "pinn_projection", 9, params);

  if (params.BARYON)
  {
    writeScalarToFile(baryonDensity, "nB", params);
    writeScalarToFileProjection(baryonDensity, "nB_projection", params);
    writeVectorToFile(baryonDiffusion, "Vt", 0, params);
    writeVectorToFile(baryonDiffusion, "Vx", 1, params);
    writeVectorToFile(baryonDiffusion, "Vy", 2, params);
    writeVectorToFile(baryonDiffusion, "Vn", 3, params);
  }

  //free the memory
  free2dArray(stressTensor, 10);
  free(energyDensity);
  free2dArray(flowVelocity, 4);
  free(pressure);
  free(bulkPressure);
  free2dArray(shearTensor, 10);

  if (params.BARYON)
  {
    free2dArray(baryonCurrent, 4);
    free(baryonDensity);
    free2dArray(baryonDiffusion, 4);
  }

  if (PRINT_SCREEN) printf("Done... Goodbye!\n");
}
