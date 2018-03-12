#pragma once
#include <unistd.h>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <string>
#include "Parameter.h"
#include <math.h>
void writeScalarToFile(double *var, char name[255], parameters params)
{
  int DIM_X = params.DIM_X;
  int DIM_Y = params.DIM_Y;
  int DIM_ETA = params.DIM_ETA;
  double DX = params.DX;
  double DY = params.DY;
  double DETA = params.DETA;
  char cwd[1024];
  getcwd(cwd, sizeof(cwd));
  std::ofstream myfile;
  char filename[255] = "";
  sprintf(filename, "output/%s.dat", name);
  myfile.open(filename);
  for (int ix = 0; ix < DIM_X; ix++)
  {
    for (int iy = 0; iy < DIM_Y; iy++)
    {
      for (int ieta = 0; ieta < DIM_ETA; ieta++)
      {
        double x = (double)ix * DX  - (((double)(DIM_X-1)) / 2.0 * DX);
        x = DX * roundf(x / DX);
        double y = (double)iy * DY  - (((double)(DIM_Y-1)) / 2.0 * DY);
        y = DY * roundf(y / DY);
        double eta = (double)ieta * DETA  - (((double)(DIM_ETA-1)) / 2.0 * DETA);
        eta = DETA * roundf(eta / DETA);

        int is = (DIM_Y * DIM_ETA) * ix + (DIM_ETA) * iy + ieta; //the column packed index spanning x, y, z

        myfile << x << " " << y << " " << eta << " " << var[is] << "\n";
      }
    }
  }
  myfile.close();
}

void writeVectorToFile(double **var, char name[255], int idx, parameters params)
{
  int DIM_X = params.DIM_X;
  int DIM_Y = params.DIM_Y;
  int DIM_ETA = params.DIM_ETA;
  double DX = params.DX;
  double DY = params.DY;
  double DETA = params.DETA;
  char cwd[1024];
  getcwd(cwd, sizeof(cwd));
  std::ofstream myfile;
  char filename[255] = "";
  sprintf(filename, "output/%s.dat", name);
  myfile.open(filename);
  for (int ix = 0; ix < DIM_X; ix++)
  {
    for (int iy = 0; iy < DIM_Y; iy++)
    {
      for (int ieta = 0; ieta < DIM_ETA; ieta++)
      {
        double x = (double)ix * DX  - (((double)(DIM_X-1)) / 2.0 * DX);
        x = DX * roundf(x / DX); //rounding for regularly spaced values
        double y = (double)iy * DY  - (((double)(DIM_Y-1)) / 2.0 * DY);
        y = DY * roundf(y / DY);
        double eta = (double)ieta * DETA  - (((double)(DIM_ETA-1)) / 2.0 * DETA);
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
void writeScalarToFileProjection(double *var, char name[255], parameters params)
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

void writeVectorToFileProjection(double **var, char name[255], int idx, parameters params)
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

void readDensityFile(double *density, char name[255], parameters params)
{
  int DIM_X = params.DIM_X;
  int DIM_Y = params.DIM_Y;
  int DIM_ETA = params.DIM_ETA;
  double DX = params.DX;
  double DY = params.DY;
  double DETA = params.DETA;
  double xmin = (-1.0) * ((double)(DIM_X-1) / 2.0) * DX;
  double ymin = (-1.0) * ((double)(DIM_Y-1) / 2.0) * DY;
  double etamin = (-1.0) * ((double)(DIM_ETA-1) / 2.0) * DETA;
  double x, y, eta, value;

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
  double dummydouble;

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
    fscanf(fileIn, "%s\t%f\n", dummyChar, &dummydouble);
    params.ETA_WIDTH = dummydouble;
    fscanf(fileIn, "%s\t%f\n", dummyChar, &dummydouble);
    params.ETA_FLAT = dummydouble;
    fscanf(fileIn, "%s\t%f\n", dummyChar, &dummydouble);
    params.SIGMA = dummydouble;
    fscanf(fileIn, "%s\t%f\n", dummyChar, &dummydouble);
    params.SIGMA_B = dummydouble;
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
    fscanf(fileIn, "%s\t%f\n", dummyChar, &dummydouble);
    params.DX = dummydouble;
    fscanf(fileIn, "%s\t%f\n", dummyChar, &dummydouble);
    params.DY = dummydouble;
    fscanf(fileIn, "%s\t%f\n", dummyChar, &dummydouble);
    params.DETA = dummydouble;
    fscanf(fileIn, "%s\t%f\n", dummyChar, &dummydouble);
    params.DRAP = dummydouble;
    fscanf(fileIn, "%s\t%f\n", dummyChar, &dummydouble);
    params.DTAU = dummydouble;
    fscanf(fileIn, "%s\t%f\n", dummyChar, &dummydouble);
    params.TAU0 = dummydouble;
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
