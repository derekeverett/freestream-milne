#pragma once
#include "Parameter.h"
#include <math.h>
void calculatePressure(float *energyDensity, float *baryonDensity, float *pressure, parameters params)
{
  int EOS_TYPE = params.EOS_TYPE;
  int DIM = params.DIM;
  //conformal eqn of state
  if (EOS_TYPE == 1)
  {
    #pragma omp parallel for
    for (int is = 0; is < DIM; is++)
    {
      pressure[is] = energyDensity[is] / 3.0;
    }
  }
  //parameterization from Wuppertal-Budapest collaboration, taken from cpu-vh/.../EquationOfState.cpp
  //requires zero baryon density
  else if (EOS_TYPE == 2)
  {
    float a0 = -0.25181736420168666;
    float a1 = 9737.845799644809;
    float a2 = 1.077580993288114e6;
    float a3 = 3.1729694865420084e6;
    float a4 = 1.6357487344679043e6;
    float a5 = 334334.4309240126;
    float a6 = 41913.439282708554;
    float a7 = 6340.448389300905;
    float a8 = 141.5073484468774;
    float a9 = 0.7158279081255019;
    float a10 = 0.0009417586777847889;
    float a11 = 3.1188455176941583e-7;
    float a12 = 1.9531729608963267e-11;

    float b0 = 45829.44617893836;
    float b1 = 4.0574329080826794e6;
    float b2 = 2.0931169138134286e7;
    float b3 = 1.3512402226067686e7;
    float b4 = 1.7851642641834426e6;
    float b5 = 278581.2989342773;
    float b6 = 26452.34905933697;
    float b7 = 499.04919730607065;
    float b8 = 2.3405487982094204;
    float b9 = 0.002962497695527404;
    float b10 = 9.601103399348206e-7;
    float b11 = 5.928138360995685e-11;
    float b12 = 3.2581066229887368e-18;

    #pragma omp parallel for
    for (int is = 0; is < DIM; is++)
    {
      float e = energyDensity[is];
      float e1 = e;
      float e2 = e*e;
      float e3 = e2*e;
      float e4 = e3*e;
      float e5 = e4*e;
      float e6 = e5*e;
      float e7 = e6*e;
      float e8 = e7*e;
      float e9 = e8*e;
      float e10 = e9*e;
      float e11 = e10*e;
      float e12 = e11*e;
      float a = (float)fmaf(a12,e12,fmaf(a11,e11,fmaf(a10,e10,fmaf(a9,e9,fmaf(a8,e8,fmaf(a7,e7,fmaf(a6,e6,fmaf(a5,e5,fmaf(a4,e4,fmaf(a3,e3,fmaf(a2,e2,fmaf(a1,e1,a0))))))))))));
      float b = (float)fmaf(b12,e12,fmaf(b11,e11,fmaf(b10,e10,fmaf(b9,e9,fmaf(b8,e8,fmaf(b7,e7,fmaf(b6,e6,fmaf(b5,e5,fmaf(b4,e4,fmaf(b3,e3,fmaf(b2,e2,fmaf(b1,e1,b0))))))))))));
      pressure[is] = a / b;
    }
  }
  //import lattice qcd tables for finite baryon density and interpolate
  else if (EOS_TYPE == 3)
  {
    #pragma omp parallel for
    for (int is = 0; is < DIM; is++)
    {
      pressure[is] = 0.0; //fix this! just a placeholder for an interpolation of tables
    }
  }
}
