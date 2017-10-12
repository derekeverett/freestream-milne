#include <math.h>
#include <stdlib.h>
#include <fstream>

#define THETA_FUNCTION(X) ((double)X < (double)0 ? (double)0 : (double)1)

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
    int ix = is / (DIM_Y * DIM_ETA);
    int iy = (is - (DIM_Y * DIM_ETA * ix))/ DIM_ETA;
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
    int ix = is / (DIM_Y * DIM_ETA);
    int iy = (is - (DIM_Y * DIM_ETA * ix))/ DIM_ETA;
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
    int ix = is / (DIM_Y * DIM_ETA);
    int iy = (is - (DIM_Y * DIM_ETA * ix))/ DIM_ETA;
    int ieta = is - (DIM_Y * DIM_ETA * ix) - (DIM_ETA * iy);

    //does it work for even number of points?
    float x = (float)ix * DX  - ((float)(DIM_X-1)) / 2.0 * DX;
    float y = (float)iy * DY  - ((float)(DIM_Y-1)) / 2.0 * DY;
    float eta = (float)ieta * DETA  - ((float)(DIM_ETA-1)) / 2.0 * DETA;

    density[is] = ((float)rand() / RAND_MAX) * exp(-(1.0 / bx) * (x * x)) * exp(-(1.0 / by) * (y * y)) * exp(-(1.0 / beta) * (eta * eta));
  }
}

void readEnergyDensitySuperMCBlock(float *density, float etaWidth, float etaFlat)
{
  //first read in the transverse profile from superMC block data format
  float temp = 0.0;
  std::ifstream superMCFile;
  superMCFile.open("initial_superMC_ed/data_12_2760_2.5fm_51pts/ed_event_7_block.dat");
  for (int ix = 0; ix < DIM_X; ix++)
  {
    for (int iy = 0; iy < DIM_Y; iy++)
    {
      superMCFile >> temp;
      for (int ieta = 0; ieta < DIM_ETA; ieta++) //copy the same value for all eta, then we will multiply by eta dependent function
      {
        int is = (DIM_Y * DIM_ETA) * ix + (DIM_ETA) * iy + ieta; //the column packed index spanning x, y, z
        density[is] = temp;
      }
    }
  }
  superMCFile.close();

  //now multiply by an eta-dependent profile; etaWidth is the width of the eta profile
  for (int is = 0; is < DIM; is++)
  {
    int ix = is / (DIM_Y * DIM_ETA);
    int iy = (is - (DIM_Y * DIM_ETA * ix))/ DIM_ETA;
    int ieta = is - (DIM_Y * DIM_ETA * ix) - (DIM_ETA * iy);

    float eta = (float)ieta * DETA  - ((float)(DIM_ETA-1)) / 2.0 * DETA;
    //here we use a the same profile as GPU-VH (see arXiv:1608.06577v1 p. 38)
    float arg = (-1.0) * (abs(eta) - etaFlat) * (abs(eta) - etaFlat) / (2.0 * etaWidth * etaWidth);
    arg = arg * THETA_FUNCTION(abs(eta) - etaFlat);
    density[is] = density[is] * exp(arg);
  }
}
