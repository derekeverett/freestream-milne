#pragma once
#include <math.h>
#include "Parameter.h"

#ifdef _OPENACC
#include <accelmath.h>
#endif

#define PI 3.141592654f

//#pragma acc routine //build a copy of function to run on device
void freeStream(float **density, float ***shiftedDensity, parameters params)
{
  int DIM_X = params.DIM_X;
  int DIM_Y = params.DIM_Y;
  int DIM_ETA = params.DIM_ETA;
  int DIM_RAP = params.DIM_RAP;
  int DIM_PHIP = params.DIM_PHIP;
  int DIM = params.DIM;
  float DX = params.DX;
  float DY = params.DY;
  float DETA = params.DETA;
  float DRAP = params.DRAP;
  float TAU0 = params.TAU0;
  float TAU = params.TAU;

  float xmin = (-1.0) * ((float)(DIM_X-1) / 2.0) * DX;
  float ymin = (-1.0) * ((float)(DIM_Y-1) / 2.0) * DY;
  float etamin = (-1.0) * ((float)(DIM_ETA-1) / 2.0) * DETA;
  float rapmin = (-1.0) * ((float)(DIM_RAP-1) / 2.0) * DRAP;

  #pragma omp parallel for simd
  //#pragma acc parallel loop
  for (int is = 0; is < DIM; is++)
  {
    int ix = is / (DIM_Y * DIM_ETA);
    int iy = (is - (DIM_Y * DIM_ETA * ix))/ DIM_ETA;
    int ieta = is - (DIM_Y * DIM_ETA * ix) - (DIM_ETA * iy);

    float x = (float)ix * DX  + xmin;
    float y = (float)iy * DY  + ymin;
    float eta = (float)ieta * DETA  + etamin;
    for (int irap = 0; irap < DIM_RAP; irap++)
    {
      float rap = (float)irap * DRAP + rapmin;

      for (int iphip = 0; iphip < DIM_PHIP; iphip++)
      {
        float phip = float(iphip) * (2.0 * PI) / float(DIM_PHIP);

        //can these trig and hypertrig functions be tabulated ahead of time?
        //check these for correctness
        float eta_new = asinh( (TAU0 / TAU) * sinh(eta - rap) ) + rap;
        //float eta_new = asinh( (TAU / TAU0) * sinh(eta - rap) ) + rap; //old formula (before typo found)
        float x_new = x - cos(phip) * (TAU * cosh(rap - eta_new) - TAU0 * cosh(rap - eta));
        float y_new = y - sin(phip) * (TAU * cosh(rap - eta_new) - TAU0 * cosh(rap - eta));

        int ix_new = (int)round((x_new - xmin) / DX);
        int iy_new = (int)round((y_new - ymin) / DY);
        int ieta_new = (int)round((eta_new - etamin) / DETA);

        int is_new = (DIM_Y * DIM_ETA * ix_new) + (DIM_ETA * iy_new) + ieta_new;

        //prevent from going out of array bounds
        //note this may be causing problems! what happens when it goes out of array bounds?
        if ((ix_new >= 0) && (ix_new < DIM_X) && (iy_new >= 0) && (iy_new < DIM_Y) && (ieta_new >= 0) && (ieta_new < DIM_ETA))
        {
          shiftedDensity[is][irap][iphip] = density[is_new][irap];
        }
      }
    }
  }
}
//this creates the initial G^(tau,tau) function, a function of spatial coordinates and rapidity
//in the special case if 2+1D, this calculates F^(tau,tau) in the notation of (PRC 91, 064906)
//rapidity dependence is determined by the assumption for rapidity dependence of the initial distribution function
void convertInitialDensity(float *initialEnergyDensity, float **density, parameters params)
{
  float SIGMA = params.SIGMA;
  int DIM = params.DIM;
  float TAU0 = params.TAU0;
  //int DIM_X = params.DIM_X;
  int DIM_Y = params.DIM_Y;
  int DIM_ETA = params.DIM_ETA;

  int DIM_RAP = params.DIM_RAP;
  float DRAP = params.DRAP;
  float DETA = params.DETA;

  float n = (sqrt(PI) / 2.0) * SIGMA * (1.0 + exp(SIGMA * SIGMA)); //the integral over cosh^2 * exp()
  float norm_factor = 1.0 / (2.0 * PI * n); //the normalization constant relating the intial energy density to the intial density profile G(tilde)^(tau,tau)

  float rapmin = (-1.0) * ((float)(DIM_RAP-1) / 2.0) * DRAP;
  float etamin = (-1.0) * ((float)(DIM_ETA-1) / 2.0) * DETA;

  if (DIM_ETA == 1) //catch the special case of 2+1D freestreaming; note DIM_RAP must also be 1 ! the normalization with SIGMA -> 0 does not generalize?
  {
    for (int is = 0; is < DIM; is++)
    {
      int ix = is / (DIM_Y);
      int iy = (is - (DIM_Y * ix));

      for (int irap = 0; irap < DIM_RAP; irap++)
      {
        //density[is][irap] = initialEnergyDensity[is] * (TAU0 / (2.0 * PI)); //this is initial F^(tau,tau) in the notation of (PRC 91, 064906)
        density[is][irap] = initialEnergyDensity[is] / (2.0 * PI); //this is initial F^(tau,tau) in the notation of (PRC 91, 064906)

      }
    }
  }

  else
  {
    for (int is = 0; is < DIM; is++)
    {
      int ix = is / (DIM_Y * DIM_ETA);
      int iy = (is - (DIM_Y * DIM_ETA * ix))/ DIM_ETA;
      int ieta = is - (DIM_Y * DIM_ETA * ix) - (DIM_ETA * iy);

      float eta = (float)ieta * DETA  + etamin;

      for (int irap = 0; irap < DIM_RAP; irap++)
      {
        float rap = (float)irap * DRAP + rapmin;
        float rap_factor = cosh(eta - rap) * cosh(eta - rap) * exp( (-1.0) * (eta - rap) * (eta - rap) / (SIGMA * SIGMA) );
        density[is][irap] = initialEnergyDensity[is] * norm_factor * rap_factor; //this is initial G^(tau,tau)
      }
    }
  }
}
//this creates the initial J^(tau) function, a function of spatial coordinates and rapidity
//rapidity dependence is determined by the assumption for rapidity dependence of the initial baryon distribution function
void convertInitialChargeDensity(float *initialChargeDensity, float **chargeDensity, parameters params)
{
  float SIGMA_B = params.SIGMA_B;
  int DIM = params.DIM;
  //int DIM_X = params.DIM_X;
  int DIM_Y = params.DIM_Y;
  int DIM_ETA = params.DIM_ETA;

  int DIM_RAP = params.DIM_RAP;
  float DRAP = params.DRAP;
  float DETA = params.DETA;

  float n = sqrt(2.0 * PI) * SIGMA_B * exp(SIGMA_B * SIGMA_B / 2.0); //the integral over cosh * exp()
  float norm_factor = 1.0 / (2.0 * PI * n); //the normalization constant relating the intial baryon density to the intial charge density profile J(tilde)^(tau)

  float rapmin = (-1.0) * ((float)(DIM_RAP-1) / 2.0) * DRAP;
  float etamin = (-1.0) * ((float)(DIM_ETA-1) / 2.0) * DETA;

  for (int is = 0; is < DIM; is++)
  {
    int ix = is / (DIM_Y * DIM_ETA);
    int iy = (is - (DIM_Y * DIM_ETA * ix))/ DIM_ETA;
    int ieta = is - (DIM_Y * DIM_ETA * ix) - (DIM_ETA * iy);

    float eta = (float)ieta * DETA  + etamin;

    for (int irap = 0; irap < DIM_RAP; irap++)
    {
      float rap = (float)irap * DRAP + rapmin;
      float rap_factor = cosh(eta - rap) * exp((-1.0) * (eta - rap) * (eta - rap) / (SIGMA_B * SIGMA_B));
      chargeDensity[is][irap] = initialChargeDensity[is] * norm_factor * rap_factor; //this is initial J^(tau)
    }

  }
}
