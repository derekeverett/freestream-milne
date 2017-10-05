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
//using namespace std;

int main(void)
{
  printf("Welcome to freestream-milne\n");
  printf("Parameters are ...\n");
  printf("(DIM_X, DIM_Y, DIM_ETA) = (%d, %d, %d)\n", DIM_X, DIM_Y, DIM_ETA);
  printf("(DX, DY, DETA, DTAU) = (%.2f, %.2f, %.2f, %.2f)\n", DX, DY, DETA, DTAU);

  //allocate and initialize memory
  printf("Allocating memory\n");
  //the initial energy density spatial profile
  float *initialEnergyDensity;
  initialEnergyDensity = (float *)calloc(DIM, sizeof(float));
  //the initial density G(tilde)^(tau,tau) at time t0
  float *density;
  density = (float *)calloc(DIM, sizeof(float));
  //the shited density profile G^(tau,tau) at time t
  float ***shiftedDensity;
  shiftedDensity = calloc3dArray(shiftedDensity, DIM, DIM_RAP, DIM_PHIP);
  //the ten independent components of the stress tensor
  float **stressTensor;
  stressTensor = calloc2dArray(stressTensor, 10, DIM);
  //a table containing 10 rows for 10 independent combinations of p_(mu)p_(nu)
  float ***hypertrigTable;
  hypertrigTable = calloc3dArray(trigTable, 11, DIM_RAP, DIM_PHIP);

  //variables to store the hydrodynamic variables after the Landau matching is performed
  //the energy density
  float *energyDensity;
  energyDensity = (float *)calloc(DIM, sizeof(float));
  //the flow velocity
  float **flowVelocity;
  flowVelocity = calloc2dArray(flowVelocity, 4, DIM);
  //the pressure
  float *pressure;
  pressure = (float *)calloc(DIM, sizeof(float));
  float *bulkPressure;
  bulkPressure = (float *)calloc(DIM, sizeof(float));
  float **shearTensor;
  shearTensor = calloc2dArray(shearTensor, 6, DIM); //calculate 6 components, can check tracelessness for accuracy

  //initialize energy density
  printf("setting initial conditions on energy density\n");
  //initializeGauss(initialEnergyDensity, 1.0);
  //initializeMCGauss(initialEnergyDensity, 1.0);
  //initializeEllipticalGauss(initialEnergyDensity, 0.5, 1.0, 3.0);
  initializeEllipticalMCGauss(initialEnergyDensity, 0.5, 1.0, 3.0);

  //write initial profile to file
  writeScalarToFile(initialEnergyDensity, "initial_e");
  writeScalarToFileProjection(initialEnergyDensity, "initial_e_projection");

  //read in the initial energy density profile (from file)
  //readInitialEnergyDensity(initialEnergyDensity);

  //convert the energy density profile into the initial density profile to be streamed - just a normalization
  convertInitialDensity(initialEnergyDensity, density);

  //perform the free streaming time-update step
  //pretabulate trig and hypertrig functions before this step to save time
  printf("performing the free streaming time step\n");
  double sec;
  sec = omp_get_wtime();
  freeStream(density, shiftedDensity);
  sec = omp_get_wtime() - sec;
  printf("Free streaming took %f seconds\n", sec);

  //Landau matching to find the components of energy-momentum tensor
  printf("Landau matching to find hydrodynamic variables\n");

  //printf("calculating trig table\n");
  sec = omp_get_wtime();
  calculateTrigTable(hypertrigTable);
  sec = omp_get_wtime() - sec;
  //printf("calculating trig table took %f seconds\n", sec);

  //calculate the ten independent components of the stress tensor by integrating over momentum angles
  printf("calculating independent components of stress tensor\n");
  sec = omp_get_wtime();
  calculateStressTensor(stressTensor, shiftedDensity, hypertrigTable);
  sec = omp_get_wtime() - sec;
  printf("calculating stress tensor took %f seconds\n", sec);

  //solve the eigenvalue problem for the energy density and flow velocity
  printf("solving eigenvalue problem for energy density and flow velocity\n");
  sec = omp_get_wtime();
  solveEigenSystem(stressTensor, energyDensity, flowVelocity);
  sec = omp_get_wtime() - sec;
  printf("solving eigenvalue problem took %f seconds\n", sec);

  calculatePressure(energyDensity, pressure);
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

  //free the memory
  free(initialEnergyDensity);
  free(density);
  free3dArray(shiftedDensity, DIM, DIM_RAP);
  free2dArray(stressTensor, 10);
  free3dArray(hypertrigTable, 11, DIM_RAP);

  free(energyDensity);
  free2dArray(flowVelocity, 4);
  free(pressure);
  free(bulkPressure);
  free2dArray(shearTensor, 6);

  printf("Done... Goodbye!\n");
}
