//trigTable is a table with 10 rows for ten combinations or p^(mu)p_(nu) normalized by the energy
#pragma once
#include <stdio.h>
//#include <math.h>
#include "Parameter.h"
#define PI 3.141592654f

void calculateBulkInvReynolds(float *pressure, float *bulkPressure, float *R_Pi_Inv, parameters params)
{
  int DIM = params.DIM;
  float TAU = params.TAU;
  #pragma omp parallel for
  for (int is = 0; is < DIM; is++)
  {
    float p;
    if (pressure[is] < 1.0e-5) p = 1.0e-5;
    else p = pressure[is];
    R_Pi_Inv[is] = bulkPressure[is] / p;
  }
}

void calculateShearInvReynolds(float *energyDensity, float *pressure, float **shearTensor, float *R_pimunu_Inv, parameters params)
{
  int DIM = params.DIM;
  float TAU = params.TAU;
  float tau2 = TAU*TAU;
  #pragma omp parallel for
  for (int is = 0; is < DIM; is++)
  {
    float num = shearTensor[0][is]*shearTensor[0][is] - 2.0 * (shearTensor[1][is]*shearTensor[1][is] + shearTensor[2][is]*shearTensor[2][is] + tau2*shearTensor[3][is]*shearTensor[3][is])
    + (shearTensor[4][is]*shearTensor[4][is] + shearTensor[7][is]*shearTensor[7][is] + tau2*tau2*shearTensor[9][is]*shearTensor[9][is])
    + 2.0 * (shearTensor[5][is]*shearTensor[5][is] + tau2*shearTensor[6][is]*shearTensor[6][is] + tau2*shearTensor[8][is]*shearTensor[8][is]);

    float den;
    if (energyDensity[is] < 1.0e-5) den = 1.0e-5;
    else den = energyDensity[is]*energyDensity[is] + 3.0*pressure[is]*pressure[is];
    R_pimunu_Inv[is] = sqrt(fabs(num)) / sqrt(fabs(den));
  }
}

void calculate_pi_dot_u(float **flowVelocity, float **shearTensor, float *pi_dot_u_tau, float *pi_dot_u_x, float *pi_dot_u_y, float *pi_dot_u_eta, parameters params)
{
  int DIM = params.DIM;
  float tau = params.TAU;
  #pragma omp parallel for
  for (int is = 0; is < DIM; is++)
  {
    //contravariant flow u^mu
    float ut = flowVelocity[0][is];
    float ux = flowVelocity[1][is];
    float uy = flowVelocity[2][is];
    float un = flowVelocity[3][is];
    
    //covariant flow u_mu
    ux = -ux;
    uy = -uy;
    un = -tau * tau * un;
    
    //contravariant shear stress
    float pitt = shearTensor[0][is];
    float pitx = shearTensor[1][is];
    float pity = shearTensor[2][is];
    float pitn = shearTensor[3][is];
    float pixx = shearTensor[4][is];
    float pixy = shearTensor[5][is];
    float pixn = shearTensor[6][is];
    float piyy = shearTensor[7][is];
    float piyn = shearTensor[8][is];
    float pinn = shearTensor[9][is];

    pi_dot_u_tau[is] = pitt*ut + pitx*ux + pity*uy + pitn*un;
    pi_dot_u_x[is] = pitx*ut + pixx*ux + pixy*uy + pixn*un;
    pi_dot_u_y[is] = pity*ut + pixy*ux + piyy*uy + piyn*un;
    pi_dot_u_eta[is] = pitn*ut + pixn*ux + piyn*uy + pinn*un;

  }
}

void calculate_pi_mu_mu(float **shearTensor, float *pi_mu_mu, parameters params)
{
  int DIM = params.DIM;
  float tau = params.TAU;
  #pragma omp parallel for
  for (int is = 0; is < DIM; is++)
  {
    
    //contravariant shear stress
    float pitt = shearTensor[0][is];
    float pixx = shearTensor[4][is];
    float piyy = shearTensor[7][is];
    float pinn = shearTensor[9][is];
    
    //mixed indices pi^mu_mu
    pixx = -pixx;
    piyy = -piyy;
    pinn = -tau * tau * pinn;
    
    pi_mu_mu[is] = pitt + pixx + piyy + pinn;

  }
}
