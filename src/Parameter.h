#pragma once
struct parameters
{
  int OUTPUTFORMAT;
  int BARYON;
  int IC_ENERGY;
  int IC_BARYON;
  float ETA_WIDTH;
  float ETA_FLAT;
  float SIGMA;
  float SIGMA_B;
  int DIM_X;
  int DIM_Y;
  int DIM_ETA;
  int DIM_RAP;
  int DIM_PHIP;
  float DX;
  float DY;
  float DETA;
  float DRAP;
  float DTAU;
  float TAU0;
  int NT;
  int EOS_TYPE;
  float E_FREEZE;
  int VISCOUS_MATCHING;
  int E_DEP_FS;
  float E_R;
  float TAU_R;
  float ALPHA;
  //these are computed based on the chosen parameters above; they are constrained
  int DIM;
  float TAU;
};
