#pragma once
#include <math.h>
#include "Parameter.h"

#ifdef _OPENACC
#include <accelmath.h>
#endif

#define PI 3.141592654f


float linearInterp3D(float x0, float x1, float x2,
                      float a000, float a100, float a010, float a001,
                      float a110, float a101, float a011, float a111)
{
  float result = 0.0;
  result = ( (1-x0) * (1-x1) * (1-x2) * a000 )
            + ( (x0) * (1-x1) * (1-x2) * a100 )
            + ( (1-x0) * (x1) * (1-x2) * a010 )
            + ( (1-x0) * (1-x1) * (x2) * a001 )
            + ( (x0) * (x1) * (1-x2) * a110 )
            + ( (x0) * (1-x1) * (x2) * a101 )
            + ( (1-x0) * (x1) * (x2) * a011 )
            + ( (x0) * (x1) * (x2)  * a111 );

  return result;
}

float linearInterp2D(float x0, float x1,
                      float a00, float a10, float a01, float a11)
{
  float result = 0.0;
  result = ( (1-x0) * (1-x1) * a00 )
            + ( (x0) * (1-x1) * a10 )
            + ( (1-x0) * (x1) * a01 )
            + ( (x0) * (x1) * a11 );

  return result;
}

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
  //float DRAP = params.DRAP;
  float TAU0 = params.TAU0;
  float TAU = params.TAU;

  float xmin = (-1.0) * ((float)(DIM_X-1) / 2.0) * DX;
  float ymin = (-1.0) * ((float)(DIM_Y-1) / 2.0) * DY;
  float etamin = (-1.0) * ((float)(DIM_ETA-1) / 2.0) * DETA;
  //float rapmin = (-1.0) * ((float)(DIM_RAP-1) / 2.0) * DRAP;

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
      //float rap = (float)irap * DRAP + rapmin;

      //try evaluating at values of rapidity y centered around y ~= eta
      //if (DIM_ETA > 1) rap = rap + eta;

      //w is an integration variable on the domain (-1,1) - careful not to include endpoints (nans)
      float w =  -.9975 + (float)irap * (1.995 / (float)(DIM_RAP - 1));
      float rap = eta + tan((PI / 2.0) * w );

      for (int iphip = 0; iphip < DIM_PHIP; iphip++)
      {
        float phip = float(iphip) * (2.0 * PI) / float(DIM_PHIP);

        //can these trig and hypertrig functions be tabulated ahead of time?
        //check these for correctness
        //float eta_new = asinh( (TAU0 / TAU) * sinh(eta - rap) ) + rap;
        float eta_new;
        if (DIM_ETA == 1) eta_new = 0.0;
        else if (DIM_ETA > 1) eta_new = asinh( (TAU / TAU0) * sinh(eta - rap) ) + rap; //old formula works
        float x_new = x - cos(phip) * (TAU * cosh(rap - eta_new) - TAU0 * cosh(rap - eta));
        float y_new = y - sin(phip) * (TAU * cosh(rap - eta_new) - TAU0 * cosh(rap - eta));

        //get fractions for interpolation routine
        //right now using a linear interpolation, which is gauranteed positive definite
        //could try cubic spline interpolation?
        float ix_new = (x_new - xmin) / DX;
        float iy_new = (y_new - ymin) / DY;
        float ieta_new = (eta_new - etamin) / DETA;

        float ix_new_f = floor(ix_new);
        float iy_new_f = floor(iy_new);
        float ieta_new_f = floor(ieta_new);

        float x_frac = ix_new - ix_new_f;
        float y_frac = iy_new - iy_new_f;
        float eta_frac = ieta_new - ieta_new_f;

        //prevent from going out of array bounds!
        if ( (ix_new_f >= 1) && (ix_new_f < DIM_X - 1) && (iy_new_f >= 1) && (iy_new_f < DIM_Y - 1) )
        {
          float interp = 0.0;

          if (DIM_ETA == 1)
          {
            int is_new_00 = (DIM_Y * DIM_ETA * (int)ix_new_f) + (DIM_ETA * (int)iy_new_f) + (int)ieta_new_f;
            int is_new_10 = (DIM_Y * DIM_ETA * ( (int)ix_new_f + 1) ) + (DIM_ETA * (int)iy_new_f) + (int)ieta_new_f;
            int is_new_01 = (DIM_Y * DIM_ETA * (int)ix_new_f) + (DIM_ETA * ( (int)iy_new_f + 1) ) + (int)ieta_new_f;
            int is_new_11 = (DIM_Y * DIM_ETA * ( (int)ix_new_f + 1) ) + (DIM_ETA * ( (int)iy_new_f + 1) ) + (int)ieta_new_f;

            float a00 = density[is_new_00][irap];
            float a10 = density[is_new_10][irap];
            float a01 = density[is_new_01][irap];
            float a11 = density[is_new_11][irap];

            interp = linearInterp2D(x_frac, y_frac, a00, a10, a01, a11);
          }

          else if ( (DIM_ETA > 1) && (ieta_new_f >= 1) && (ieta_new_f < DIM_ETA - 1) )
          {
            int is_new_000 = (DIM_Y * DIM_ETA * (int)ix_new_f) + (DIM_ETA * (int)iy_new_f) + (int)ieta_new_f;
            int is_new_100 = (DIM_Y * DIM_ETA * ( (int)ix_new_f + 1) ) + (DIM_ETA * (int)iy_new_f) + (int)ieta_new_f;
            int is_new_010 = (DIM_Y * DIM_ETA * (int)ix_new_f) + (DIM_ETA * ( (int)iy_new_f + 1) ) + (int)ieta_new_f;
            int is_new_110 = (DIM_Y * DIM_ETA * ( (int)ix_new_f + 1) ) + (DIM_ETA * ( (int)iy_new_f + 1) ) + (int)ieta_new_f;
            int is_new_001 = (DIM_Y * DIM_ETA * (int)ix_new_f) + (DIM_ETA * (int)iy_new_f) + ( (int)ieta_new_f + 1 );
            int is_new_101 = (DIM_Y * DIM_ETA * ( (int)ix_new_f + 1) ) + (DIM_ETA * (int)iy_new_f) + ( (int)ieta_new_f + 1 );
            int is_new_011 = (DIM_Y * DIM_ETA * (int)ix_new_f) + (DIM_ETA * ( (int)iy_new_f + 1) ) + ( (int)ieta_new_f + 1 );
            int is_new_111 = (DIM_Y * DIM_ETA * ( (int)ix_new_f + 1) ) + (DIM_ETA * ( (int)iy_new_f + 1) ) + ( (int)ieta_new_f + 1 );

            float a000 = density[is_new_000][irap];
            float a100 = density[is_new_100][irap];
            float a010 = density[is_new_010][irap];
            float a110 = density[is_new_110][irap];
            float a001 = density[is_new_001][irap];
            float a101 = density[is_new_101][irap];
            float a011 = density[is_new_011][irap];
            float a111 = density[is_new_111][irap];

            interp = linearInterp3D(x_frac, y_frac, eta_frac,
                                  a000, a100, a010, a001,
                                  a110, a101, a011, a111);
          }

          shiftedDensity[is][irap][iphip] = interp;
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
  //float DRAP = params.DRAP;
  float DETA = params.DETA;

  float n = sqrt(PI / 2.0) * SIGMA * (1.0 + exp(2.0 * SIGMA * SIGMA)); //the integral over cosh^2 * exp()
  float norm_factor = 1.0 / (2.0 * PI * n); //the normalization constant relating the intial energy density to the intial density profile G(tilde)^(tau,tau)

  //float rapmin = (-1.0) * ((float)(DIM_RAP-1) / 2.0) * DRAP;
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
        //float rap = (float)irap * DRAP + rapmin;

        //try evaluating at values of rapidity y centered around y ~= eta
        //rap = rap + eta;

        //w is an integration variable on the domain (-1,1) - careful not to include endpoints (nans)
        float w =  -.9975 + (float)irap * (1.995 / (float)(DIM_RAP - 1));
        float rap = eta + tan((PI / 2.0) * w );

        float rap_factor = cosh(eta - rap) * cosh(eta - rap) * exp( (-1.0) * (eta - rap) * (eta - rap) / (2.0 * SIGMA * SIGMA) );
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
