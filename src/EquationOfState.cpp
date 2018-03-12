#pragma once
#include "Parameter.h"
#include <math.h>
void calculatePressure(double *energyDensity, double *baryonDensity, double *pressure, parameters params)
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
    double a0 = -0.25181736420168666;
    double a1 = 9737.845799644809;
    double a2 = 1.077580993288114e6;
    double a3 = 3.1729694865420084e6;
    double a4 = 1.6357487344679043e6;
    double a5 = 334334.4309240126;
    double a6 = 41913.439282708554;
    double a7 = 6340.448389300905;
    double a8 = 141.5073484468774;
    double a9 = 0.7158279081255019;
    double a10 = 0.0009417586777847889;
    double a11 = 3.1188455176941583e-7;
    double a12 = 1.9531729608963267e-11;

    double b0 = 45829.44617893836;
    double b1 = 4.0574329080826794e6;
    double b2 = 2.0931169138134286e7;
    double b3 = 1.3512402226067686e7;
    double b4 = 1.7851642641834426e6;
    double b5 = 278581.2989342773;
    double b6 = 26452.34905933697;
    double b7 = 499.04919730607065;
    double b8 = 2.3405487982094204;
    double b9 = 0.002962497695527404;
    double b10 = 9.601103399348206e-7;
    double b11 = 5.928138360995685e-11;
    double b12 = 3.2581066229887368e-18;

    #pragma omp parallel for
    for (int is = 0; is < DIM; is++)
    {
      double e = energyDensity[is];
      double e1 = e;
      double e2 = e*e;
      double e3 = e2*e;
      double e4 = e3*e;
      double e5 = e4*e;
      double e6 = e5*e;
      double e7 = e6*e;
      double e8 = e7*e;
      double e9 = e8*e;
      double e10 = e9*e;
      double e11 = e10*e;
      double e12 = e11*e;
      double a = (double)fma(a12,e12,fma(a11,e11,fma(a10,e10,fma(a9,e9,fma(a8,e8,fma(a7,e7,fma(a6,e6,fma(a5,e5,fma(a4,e4,fma(a3,e3,fma(a2,e2,fma(a1,e1,a0))))))))))));
      double b = (double)fma(b12,e12,fma(b11,e11,fma(b10,e10,fma(b9,e9,fma(b8,e8,fma(b7,e7,fma(b6,e6,fma(b5,e5,fma(b4,e4,fma(b3,e3,fma(b2,e2,fma(b1,e1,b0))))))))))));
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
