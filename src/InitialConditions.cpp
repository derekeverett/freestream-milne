#include <math.h>
#include <stdlib.h>

void initializeZero(float *density)
{
  for (int is = 0; is < DIM; is++)
  {
    density[is] = 0.0;
  }
}

void initializeGauss(float *density, float b) // b is the variance ('spherically' symmetric)
{
  for (int is = 0; is < DIM; is++)
  {
    int ix = is / (DIM_Y * DIM_ETA);
    int iy = (is - (DIM_Y * DIM_ETA * ix))/ DIM_ETA;
    int ieta = is - (DIM_Y * DIM_ETA * ix) - (DIM_ETA * iy);

    //does it work for even number of points?
    float x = (float)ix * DX  - ((float)(DIM_X-1)) / 2.0 * DX;
    float y = (float)iy * DY  - ((float)(DIM_Y-1)) / 2.0 * DY;
    float eta = (float)ieta * DETA  - ((float)(DIM_ETA-1)) / 2.0 * DETA;

    density[is] = exp(-(1.0 / b) * ((x * x) + (y * y) + (eta * eta)));
  }
}

void initializeEllipticalGauss(float *density, float bx, float by, float beta) // bx is the x variance etc...
{
  for (int is = 0; is < DIM; is++)
  {
    int ix = is / (DIM_Y * DIM_Z);
    int iy = (is - (DIM_Y * DIM_Z * ix))/ DIM_Z;
    int ieta = is - (DIM_Y * DIM_ETA * ix) - (DIM_ETA * iy);

    //does it work for even number of points?
    float x = (float)ix * DX  - ((float)(DIM_X-1)) / 2.0 * DX;
    float y = (float)iy * DY  - ((float)(DIM_Y-1)) / 2.0 * DY;
    float eta = (float)ieta * DETA  - ((float)(DIM_ETA-1)) / 2.0 * DETA;

    density[is] = exp(-(1.0 / bx) * (x * x)) * exp(-(1.0 / by) * (y * y)) * exp(-(1.0 / beta) * (eta * eta));
  }
}

void initializeMCGauss(float * density, float b)
{
  for (int is = 0; is < DIM; is++)
  {
    int ix = is / (DIM_Y * DIM_Z);
    int iy = (is - (DIM_Y * DIM_Z * ix))/ DIM_Z;
    int ieta = is - (DIM_Y * DIM_ETA * ix) - (DIM_ETA * iy);

    //does it work for even number of points?
    float x = (float)ix * DX  - ((float)(DIM_X-1)) / 2.0 * DX;
    float y = (float)iy * DY  - ((float)(DIM_Y-1)) / 2.0 * DY;
    float eta = (float)ieta * DETA  - ((float)(DIM_ETA-1)) / 2.0 * DETA;

    density[is] = ((float)rand() / RAND_MAX) * exp(-(1.0 / b) * ((x * x) + (y * y) + (eta * eta)));
  }
}

void initializeEllipticalMCGauss(float *density, float bx, float by, float beta) // bx is the x variance etc...
{
  for (int is = 0; is < DIM; is++)
  {
    int ix = is / (DIM_Y * DIM_Z);
    int iy = (is - (DIM_Y * DIM_Z * ix))/ DIM_Z;
    int ieta = is - (DIM_Y * DIM_ETA * ix) - (DIM_ETA * iy);

    //does it work for even number of points?
    float x = (float)ix * DX  - ((float)(DIM_X-1)) / 2.0 * DX;
    float y = (float)iy * DY  - ((float)(DIM_Y-1)) / 2.0 * DY;
    float eta = (float)ieta * DETA  - ((float)(DIM_ETA-1)) / 2.0 * DETA;

    density[is] = ((float)rand() / RAND_MAX) * exp(-(1.0 / bx) * (x * x)) * exp(-(1.0 / by) * (y * y)) * exp(-(1.0 / beta) * (eta * eta));
  }
}
