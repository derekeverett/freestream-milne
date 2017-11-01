// freestream-milne

#include "Parameters.h"
#include "FreeStream.cpp"
#include "InitialConditions.cpp"
#include "LandauMatch.cpp"
#include "EquationOfState.cpp"
#include "Memory.cpp"
#include "FileIO.cpp"
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <omp.h>

int main(void)
{
  printf("Welcome to freestream-milne\n");

  //read in chosen parameters from parameters.dat. Not yet functional.
  /*

  */
  //readInParameters(INITCOND,ETA_WIDTH,ETA_FLAT,SIGMA,PI,DIM_X,DIM_Y,DIM_ETA,DIM_RAP,DIM_PHIP,DX,DY,DETA,DRAP,DTAU,TAU0,EOS_TYPE)

  printf("Parameters are ...\n");
  printf("(DIM_X, DIM_Y, DIM_ETA) = (%d, %d, %d)\n", DIM_X, DIM_Y, DIM_ETA);
  printf("(DX, DY, DETA, DTAU) = (%.2f, %.2f, %.2f, %.2f)\n", DX, DY, DETA, DTAU);
  if (EOS_TYPE == 1) printf("Using EoS : Conformal \n");
  else if (EOS_TYPE == 2) printf("Using EoS : Wuppertal-Budhapest \n");
  else if (EOS_TYPE == 3) printf("Using EoS : Lattice QCD + HRG matched \n");

  //allocate and initialize memory
  printf("Allocating memory\n");
  //the initial energy density spatial profile
  float *initialEnergyDensity;
  initialEnergyDensity = (float *)calloc(DIM, sizeof(float));

  //the initial baryon density spatial profile
  float *initialChargeDensity;
  if(BARYON) initialChargeDensity = (float *)calloc(DIM, sizeof(float));

  //the initial density G(tilde)^(tau,tau) at time tau_0
  float **density;
  density = calloc2dArray(density, DIM, DIM_RAP); // function of x,y,eta and rapidity

  //the initial density J(tilde)^(tau) at time tau_0
  float **chargeDensity;
  if(BARYON) chargeDensity = calloc2dArray(density, DIM, DIM_RAP); // function of x,y,eta and rapidity

  //initialize energy density
  printf("setting initial conditions on energy density : ");
  if (IC_ENERGY == 1)
  {
    initializeEllipticalGauss(initialEnergyDensity, 1.0, 1.0, 1.0);
    printf("Smooth Oblate Gaussian \n");
  }
  else if (IC_ENERGY == 2)
  {
    initializeEllipticalMCGauss(initialEnergyDensity, 1.0, 1.0, 1.0);
    printf("Fluctuating Oblate Gaussian \n");
  }
  else if (IC_ENERGY == 3)
  {
    readDensityFile(initialEnergyDensity, "initial_profiles/e");
    printf("Reading from energy density file in initial_profiles/ \n");
  }
  else if (IC_ENERGY == 4)
  {
    readEnergyDensitySuperMCBlock(initialEnergyDensity, ETA_WIDTH, ETA_FLAT);
    printf("Reading from energy density file in initial_profiles/ \n");
  }
  else
  {
    printf("Not a valid initial Condition... Goodbye\n");
    return 0;
  }

  if (BARYON)
  {
    //initialize baryon density in cartesian coordinates
    printf("setting initial conditions on baryon density : ");
    if (IC_BARYON == 1)
    {
      initializeEllipticalGauss(initialChargeDensity, 1.0, 1.0, 1.0);
      printf("Smooth Oblate Gaussian \n");
    }
    else if (IC_BARYON == 2)
    {
      initializeEllipticalMCGauss(initialChargeDensity, 1.0, 1.0, 1.0);
      printf("Fluctuating Oblate Gaussian \n");
    }
    else if (IC_BARYON == 3)
    {
      readDensityFile(initialChargeDensity, "initial_profiles/nB");
      printf("Reading from baryon density file in initial_profiles/ \n");
    }
    else
    {
      printf("Not a valid initial Condition... Goodbye\n");
      return 0;
    }
  }

  //write initial energy density and baryon density to file
  writeScalarToFile(initialEnergyDensity, "initial_e");
  if (BARYON) writeScalarToFile(initialChargeDensity, "initial_nB");
  writeScalarToFileProjection(initialEnergyDensity, "initial_e_projection");
  if (BARYON) writeScalarToFileProjection(initialChargeDensity, "initial_nB_projection");

  //convert the energy density profile into the initial density profile to be streamed and free memory
  convertInitialDensity(initialEnergyDensity, density);
  free(initialEnergyDensity);
  //convert the baryon density profile into the initial baryon density profile to be streamed and free memory
  if (BARYON) convertInitialChargeDensity(initialChargeDensity, chargeDensity);
  if (BARYON) free(initialChargeDensity);

  //the shifted energy density profile G^(tau,tau) at time tau
  float ***shiftedDensity;
  shiftedDensity = calloc3dArray(shiftedDensity, DIM, DIM_RAP, DIM_PHIP);

  //the shifted baryon density profile J^(tau) at time tau
  float ***shiftedChargeDensity;
  if(BARYON) shiftedChargeDensity = calloc3dArray(shiftedChargeDensity, DIM, DIM_RAP, DIM_PHIP);

  //perform the free streaming time-update step and free up memory
  //pretabulate trig and hypertrig functions before this step to save time?
  printf("performing the free streaming\n");
  double sec;
  sec = omp_get_wtime();
  freeStream(density, shiftedDensity);
  free2dArray(density, DIM);
  if (BARYON) freeStream(chargeDensity, shiftedChargeDensity);
  if (BARYON) free2dArray(chargeDensity, DIM);

  sec = omp_get_wtime() - sec;
  printf("Free streaming took %f seconds\n", sec);

  //Landau matching to find the components of energy-momentum tensor
  printf("Landau matching to find hydrodynamic variables\n");

  //the ten independent components of the stress tensor
  float **stressTensor;
  stressTensor = calloc2dArray(stressTensor, 10, DIM);

  //the four independent components of baryon current four-vector
  float **baryonCurrent;
  if(BARYON) baryonCurrent = calloc2dArray(baryonCurrent, 4, DIM);

  //a table containing 10 rows for 10 independent combinations of p_(mu)p_(nu)
  float ****hypertrigTable;
  hypertrigTable = calloc4dArray(hypertrigTable, 10, DIM_RAP, DIM_PHIP, DIM_ETA); //depends on eta because we have function of eta - y

  printf("calculating hypertrig table\n");
  sec = omp_get_wtime();
  calculateHypertrigTable(hypertrigTable);
  sec = omp_get_wtime() - sec;
  printf("calculating trig table took %f seconds\n", sec);

  //calculate the ten independent components of the stress tensor by integrating over rapidity and phi_p
  printf("calculating independent components of stress tensor\n");
  sec = omp_get_wtime();
  calculateStressTensor(stressTensor, shiftedDensity, hypertrigTable);
  free3dArray(shiftedDensity, DIM, DIM_RAP);
  sec = omp_get_wtime() - sec;
  printf("calculating stress tensor took %f seconds\n", sec);

  if (BARYON)
  {
    //calculate the four independent components of the baryon current by integrating over rapidity and phi_p
    printf("calculating independent components of baryon current\n");
    sec = omp_get_wtime();
    calculateBaryonCurrent(baryonCurrent, shiftedChargeDensity, hypertrigTable);
    free3dArray(shiftedChargeDensity, DIM, DIM_RAP);
    sec = omp_get_wtime() - sec;
    printf("calculating baryon current took %f seconds\n", sec);
  }

  //done with hypertrig table as well
  free4dArray(hypertrigTable, 10, DIM_RAP, DIM_PHIP);

  //variables to store the hydrodynamic variables after the Landau matching is performed
  //the energy density
  float *energyDensity;
  energyDensity = (float *)calloc(DIM, sizeof(float));

  //the baryon density
  float *baryonDensity;
  if(BARYON) baryonDensity = (float *)calloc(DIM, sizeof(float));

  //the flow velocity
  float **flowVelocity;
  flowVelocity = calloc2dArray(flowVelocity, 4, DIM);

  //the pressure
  float *pressure;
  pressure = (float *)calloc(DIM, sizeof(float));

  //the bulk pressure Pi
  float *bulkPressure;
  bulkPressure = (float *)calloc(DIM, sizeof(float));

  //the shear stress tensor
  float **shearTensor;
  shearTensor = calloc2dArray(shearTensor, 10, DIM); //calculate 10 components, can check tracelessness/orthogonality for accuracy

  //the baryon diffusion current vector
  float **baryonDiffusion;
  if(BARYON) baryonDiffusion = calloc2dArray(baryonDiffusion, 4, DIM);

  //solve the eigenvalue problem for the energy density and flow velocity
  printf("solving eigenvalue problem for energy density and flow velocity\n");
  sec = omp_get_wtime();
  solveEigenSystem(stressTensor, energyDensity, flowVelocity);
  sec = omp_get_wtime() - sec;
  printf("solving eigenvalue problem took %f seconds\n", sec);

  if (BARYON)
  {
    //calculate baryon density and diffusion current
    printf("calculating baryon density and diffusion current \n");
    sec = omp_get_wtime();
    calculateBaryonDensity(baryonDensity, baryonCurrent, flowVelocity);
    calculateBaryonDiffusion(baryonDiffusion, baryonCurrent, baryonDensity, flowVelocity);
    sec = omp_get_wtime() - sec;
    printf("calculating baryon density and diffusion current took %f seconds\n", sec);
  }

  calculatePressure(energyDensity, baryonDensity, pressure);
  calculateBulkPressure(stressTensor, energyDensity, pressure, bulkPressure);
  calculateShearViscTensor(stressTensor, energyDensity, flowVelocity, pressure, bulkPressure, shearTensor);

  printf("writing hydro variables to file\n");
  writeScalarToFile(energyDensity, "e");
  writeScalarToFile(pressure, "p");
  writeScalarToFile(bulkPressure, "bulk_PI");
  writeScalarToFileProjection(energyDensity, "e_projection");
  writeVectorToFile(flowVelocity, "u_x", 1);
  writeVectorToFile(flowVelocity, "u_y", 2);
  writeVectorToFile(flowVelocity, "u_eta", 3);
  writeVectorToFileProjection(flowVelocity, "u_x_projection", 1);
  writeVectorToFileProjection(flowVelocity, "u_y_projection", 2);
  writeVectorToFileProjection(flowVelocity, "u_eta_projection", 3);

  if (BARYON)
  {
    writeScalarToFile(baryonDensity, "nB");
    writeScalarToFileProjection(baryonDensity, "nB_projection");
    writeVectorToFile(baryonDiffusion, "V_x", 1);
    writeVectorToFile(baryonDiffusion, "V_y", 2);
    writeVectorToFile(baryonDiffusion, "V_eta", 3);
  }

  //free the memory
  free2dArray(stressTensor, 10);
  free(energyDensity);
  free2dArray(flowVelocity, 4);
  free(pressure);
  free(bulkPressure);
  free2dArray(shearTensor, 10);

  if (BARYON)
  {
    free2dArray(baryonCurrent, 4);
    free(baryonDensity);
    free2dArray(baryonDiffusion, 4);
  }

  printf("Done... Goodbye!\n");
}
