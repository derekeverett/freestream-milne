#pragma once
struct parameters
{
  //sets defualt values if can't find parameters.dat file
  int BARYON = 0;
  int IC_ENERGY = 1;
  int IC_BARYON = 1;
  float ETA_WIDTH = 0.5;
  float ETA_FLAT = 0.5;
  float SIGMA = 1.0;
  float SIGMA_B = 1.0;
  int DIM_X = 51;
  int DIM_Y = 51;
  int DIM_ETA = 51;
  int DIM_RAP = 51;
  int DIM_PHIP = 31;
  float DX = 0.1;
  float DY = 0.1;
  float DETA = 0.1;
  float DRAP = 0.2;
  float DTAU = 0.5;
  float TAU0 = 0.1;
  int EOS_TYPE = 1;

  int DIM = DIM_X * DIM_Y * DIM_ETA;
  float TAU = TAU0 + DTAU;
};
