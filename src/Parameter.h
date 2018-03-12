#pragma once
struct parameters
{
  int OUTPUTFORMAT;
  int BARYON;
  int IC_ENERGY;
  int IC_BARYON;
  double ETA_WIDTH;
  double ETA_FLAT;
  double SIGMA;
  double SIGMA_B;
  int DIM_X;
  int DIM_Y;
  int DIM_ETA;
  int DIM_RAP;
  int DIM_PHIP;
  double DX;
  double DY;
  double DETA;
  double DRAP;
  double DTAU;
  double TAU0;
  int EOS_TYPE;

  //these are computed based on the chosen parameters above; they are constrained
  int DIM;
  double TAU;
};
