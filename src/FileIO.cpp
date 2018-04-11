#pragma once
#include <unistd.h>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <string>
#include "Parameter.h"
#include <math.h>

#define MINX 0
#define MINY 0
#define MINETA 0
void writeScalarToFile(float *var, char name[255], parameters params)
{
  int DIM_X = params.DIM_X;
  int DIM_Y = params.DIM_Y;
  int DIM_ETA = params.DIM_ETA;
  float DX = params.DX;
  float DY = params.DY;
  float DETA = params.DETA;
  char cwd[1024];
  getcwd(cwd, sizeof(cwd));
  std::ofstream myfile;
  char filename[255] = "";
  sprintf(filename, "output/%s.dat", name);
  myfile.open(filename);
  for (int ix = MINX; ix < DIM_X - MINX; ix++)
  {
    for (int iy = MINY; iy < DIM_Y - MINY; iy++)
    {
      for (int ieta = MINETA; ieta < DIM_ETA - MINETA; ieta++)
      {
        float x = (float)ix * DX  - (((float)(DIM_X-1)) / 2.0 * DX);
        x = DX * roundf(x / DX);
        float y = (float)iy * DY  - (((float)(DIM_Y-1)) / 2.0 * DY);
        y = DY * roundf(y / DY);
        float eta = (float)ieta * DETA  - (((float)(DIM_ETA-1)) / 2.0 * DETA);
        eta = DETA * roundf(eta / DETA);

        int is = (DIM_Y * DIM_ETA) * ix + (DIM_ETA) * iy + ieta; //the column packed index spanning x, y, z

        myfile << x << " " << y << " " << eta << " " << var[is] << "\n";
      }
    }
  }
  myfile.close();
}

void writeVectorToFile(float **var, char name[255], int idx, parameters params)
{
  int DIM_X = params.DIM_X;
  int DIM_Y = params.DIM_Y;
  int DIM_ETA = params.DIM_ETA;
  float DX = params.DX;
  float DY = params.DY;
  float DETA = params.DETA;
  char cwd[1024];
  getcwd(cwd, sizeof(cwd));
  std::ofstream myfile;
  char filename[255] = "";
  sprintf(filename, "output/%s.dat", name);
  myfile.open(filename);
  for (int ix = MINX; ix < DIM_X - MINX; ix++)
  {
    for (int iy = MINY; iy < DIM_Y - MINY; iy++)
    {
      for (int ieta = MINETA; ieta < DIM_ETA - MINETA; ieta++)
      {
        float x = (float)ix * DX  - (((float)(DIM_X-1)) / 2.0 * DX);
        x = DX * roundf(x / DX); //rounding for regularly spaced values
        float y = (float)iy * DY  - (((float)(DIM_Y-1)) / 2.0 * DY);
        y = DY * roundf(y / DY);
        float eta = (float)ieta * DETA  - (((float)(DIM_ETA-1)) / 2.0 * DETA);
        eta = DETA * roundf(eta / DETA);

        int is = (DIM_Y * DIM_ETA) * ix + (DIM_ETA) * iy + ieta; //the column packed index spanning x, y, z

        myfile << x << " " << y << " " << eta << " " << var[idx][is] << "\n";
      }
    }
  }
  myfile.close();
}

//this function writes the transverse density of a variable at z = 0
// as regularly spaced values
void writeScalarToFileProjection(float *var, char name[255], parameters params)
{
  int DIM_X = params.DIM_X;
  int DIM_Y = params.DIM_Y;
  int DIM_ETA = params.DIM_ETA;
  char cwd[1024];
  getcwd(cwd, sizeof(cwd));
  std::ofstream myfile;
  char filename[255] = "";
  sprintf(filename, "output/%s.dat", name);
  myfile.open(filename);
  for (int iy = 0; iy < DIM_Y; iy++)
  {
    for (int ix = 0; ix < DIM_X; ix++)
    {
      int ieta = (DIM_ETA - 1) / 2; // at eta = 0
      int is = (DIM_Y * DIM_ETA) * ix + (DIM_ETA) * iy + ieta; //the column packed index spanning x, y, eta
      myfile << var[is] << " "; //different columns for x values
    }
    myfile << "\n"; // different rows correspond to different y values
  }
  myfile.close();
}

void writeVectorToFileProjection(float **var, char name[255], int idx, parameters params)
{
  int DIM_X = params.DIM_X;
  int DIM_Y = params.DIM_Y;
  int DIM_ETA = params.DIM_ETA;
  char cwd[1024];
  getcwd(cwd, sizeof(cwd));
  std::ofstream myfile;
  char filename[255] = "";
  sprintf(filename, "output/%s.dat", name);
  myfile.open(filename);
  for (int iy = 0; iy < DIM_Y; iy++)
  {
    for (int ix = 0; ix < DIM_X; ix++)
    {
      int ieta = (DIM_ETA - 1) / 2; //at eta = 0
      int is = (DIM_Y * DIM_ETA) * ix + (DIM_ETA) * iy + ieta; //the column packed index spanning x, y, z
      myfile << var[idx][is] << " "; //different columns for x values
    }
    myfile << "\n"; // different rows correspond to different y values
  }
  myfile.close();
}

void readDensityFile(float *density, char name[255], parameters params)
{
  int DIM_X = params.DIM_X;
  int DIM_Y = params.DIM_Y;
  int DIM_ETA = params.DIM_ETA;
  float DX = params.DX;
  float DY = params.DY;
  float DETA = params.DETA;
  float xmin = (-1.0) * ((float)(DIM_X-1) / 2.0) * DX;
  float ymin = (-1.0) * ((float)(DIM_Y-1) / 2.0) * DY;
  float etamin = (-1.0) * ((float)(DIM_ETA-1) / 2.0) * DETA;
  float x, y, eta, value;

  char filename[255] = "";
  sprintf(filename, "%s.dat", name);
  std::ifstream infile;
  infile.open(filename);
  if (!infile)
  {
    printf("Couldn't open initial profile!\n");
    exit(1);
  }
  while (infile >> x >> y >> eta >> value)
  {
    int ix = (int)round((x - xmin) / DX);
    int iy = (int)round((y - ymin) / DY);
    int ieta = (int)round((eta - etamin) / DETA);
    int is = (DIM_Y * DIM_ETA * ix) + (DIM_ETA * iy) + ieta;
    density[is] = value;
  }
  infile.close();
}

void readInParameters(struct parameters &params)
{
  char dummyChar[255];
  int dummyInt;
  //long dummyLong;
  float dummyFloat;

  FILE *fileIn;
  std::stringstream paramsStream;
  paramsStream << "freestream_input";
  fileIn = fopen(paramsStream.str().c_str(),"r");

  if (fileIn == NULL)
  {
    printf("Couldn't open parameters.dat . Using default values!\n");
  }

  else
  {
    fscanf(fileIn, "%s\t%d\n", dummyChar, &dummyInt);
    params.OUTPUTFORMAT = dummyInt;
    fscanf(fileIn, "%s\t%d\n", dummyChar, &dummyInt);
    params.BARYON = dummyInt;
    fscanf(fileIn, "%s\t%d\n", dummyChar, &dummyInt);
    params.IC_ENERGY = dummyInt;
    fscanf(fileIn, "%s\t%d\n", dummyChar, &dummyInt);
    params.IC_BARYON = dummyInt;
    fscanf(fileIn, "%s\t%f\n", dummyChar, &dummyFloat);
    params.ETA_WIDTH = dummyFloat;
    fscanf(fileIn, "%s\t%f\n", dummyChar, &dummyFloat);
    params.ETA_FLAT = dummyFloat;
    fscanf(fileIn, "%s\t%f\n", dummyChar, &dummyFloat);
    params.SIGMA = dummyFloat;
    fscanf(fileIn, "%s\t%f\n", dummyChar, &dummyFloat);
    params.SIGMA_B = dummyFloat;
    fscanf(fileIn, "%s\t%d\n", dummyChar, &dummyInt);
    params.DIM_X = dummyInt;
    fscanf(fileIn, "%s\t%d\n", dummyChar, &dummyInt);
    params.DIM_Y = dummyInt;
    fscanf(fileIn, "%s\t%d\n", dummyChar, &dummyInt);
    params.DIM_ETA = dummyInt;
    fscanf(fileIn, "%s\t%d\n", dummyChar, &dummyInt);
    params.DIM_RAP = dummyInt;
    fscanf(fileIn, "%s\t%d\n", dummyChar, &dummyInt);
    params.DIM_PHIP = dummyInt;
    fscanf(fileIn, "%s\t%f\n", dummyChar, &dummyFloat);
    params.DX = dummyFloat;
    fscanf(fileIn, "%s\t%f\n", dummyChar, &dummyFloat);
    params.DY = dummyFloat;
    fscanf(fileIn, "%s\t%f\n", dummyChar, &dummyFloat);
    params.DETA = dummyFloat;
    fscanf(fileIn, "%s\t%f\n", dummyChar, &dummyFloat);
    params.DRAP = dummyFloat;
    fscanf(fileIn, "%s\t%f\n", dummyChar, &dummyFloat);
    params.DTAU = dummyFloat;
    fscanf(fileIn, "%s\t%f\n", dummyChar, &dummyFloat);
    params.TAU0 = dummyFloat;
    fscanf(fileIn, "%s\t%d\n", dummyChar, &dummyInt);
    params.EOS_TYPE = dummyInt;

    fclose(fileIn);
  }

}
/*
void readEoSTable()
{

}
*/
