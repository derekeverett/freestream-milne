#pragma once
#include <unistd.h>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <string>
#include "Parameter.h"
#include "FSConfig.h"
#include <math.h>
#include "EquationOfState.cpp"

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
  for (int ix = 0; ix < DIM_X; ix++)
  {
    for (int iy = 0; iy < DIM_Y; iy++)
    {
      for (int ieta = 0; ieta < DIM_ETA; ieta++)
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
  for (int ix = 0; ix < DIM_X; ix++)
  {
    for (int iy = 0; iy < DIM_Y; iy++)
    {
      for (int ieta = 0; ieta < DIM_ETA; ieta++)
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
    params.NT = dummyInt;
    fscanf(fileIn, "%s\t%d\n", dummyChar, &dummyInt);
    params.EOS_TYPE = dummyInt;
    fscanf(fileIn, "%s\t%f\n", dummyChar, &dummyFloat);
    params.E_FREEZE = dummyFloat;
    fscanf(fileIn, "%s\t%d\n", dummyChar, &dummyInt);
    params.VISCOUS_MATCHING = dummyInt;

    fscanf(fileIn, "%s\t%d\n", dummyChar, &dummyInt);
    params.E_DEP_FS = dummyInt;

    fscanf(fileIn, "%s\t%f\n", dummyChar, &dummyFloat);
    params.E_R = dummyFloat;

    fscanf(fileIn, "%s\t%f\n", dummyChar, &dummyFloat);
    params.TAU_R = dummyFloat;

    fscanf(fileIn, "%s\t%f\n", dummyChar, &dummyFloat);
    params.ALPHA = dummyFloat;

    fclose(fileIn);
  }

}


// This function outputs freestreaming evolution file in binary format
void outputEvolutionDataXYEta_chun(float *energyDensity, float **flowVelocity, int itau, float tau, float tau_step, parameters params) 
{
    const float hbarc = 0.197326938;
    const float eC = params.E_FREEZE;
    const float tau0 = params.TAU0;
    
    const std::string out_name_xyeta = "evolution_all_xyeta_fs.dat";
    std::string out_open_mode;
    FILE *out_file_xyeta;
    // If it's the first timestep, overwrite the previous file
    if (tau == tau0) out_open_mode = "wb";
    else out_open_mode = "ab";
    out_file_xyeta = fopen(out_name_xyeta.c_str(), out_open_mode.c_str());
    // write out header
    const int nx        = params.DIM_X;
    const int ny        = params.DIM_Y;
    const int neta      = params.DIM_ETA;
    const float dx     = params.DX;
    const float dy     = params.DY;
    const float deta   = params.DETA;
    float xmin = (-1.0) * ((float)(nx-1) / 2.0) * dx;
    float ymin = (-1.0) * ((float)(ny-1) / 2.0) * dy;
    float etamin = (-1.0) * ((float)(neta-1) / 2.0) * deta;
    if (tau == tau0) 
    {
        const int nVar_per_cell = 10;
        float header[] = {
            static_cast<float>(tau0), 
            static_cast<float>(tau_step),
            static_cast<float>(nx), 
            static_cast<float>(dx),
            static_cast<float>(xmin),
            static_cast<float>(ny), 
            static_cast<float>(dy),
            static_cast<float>(ymin),
            static_cast<float>(neta), 
            static_cast<float>(deta),
            static_cast<float>(etamin),
            static_cast<float>(nVar_per_cell)
        };
        fwrite(header, sizeof(float), 16, out_file_xyeta);
    }
    for (int ieta = 0; ieta < neta; ieta += 1) 
    {
        for (int iy = 0; iy < ny; iy += 1) 
        {
            for (int ix = 0; ix < nx; ix += 1) 
            {
                int is = (ny * neta) * ix + (neta) * iy + ieta; //the column packed index spanning x, y, z
                float e_local    = energyDensity[is];  // 1/fm^4
                float p_local    = e_local / 3.; // (conformal EoS)
                float ux   = flowVelocity[1][is];
                float uy   = flowVelocity[2][is];
                float ueta = flowVelocity[3][is];
                // T_local is in 1/fm
                float T_local = temperatureFromEnergyDensity(e_local); // (conformal EoS)
                if (e_local*hbarc < params.E_FREEZE) continue;
                // only ouput fluid cells that are above cut-off energy density
                float ideal[] = {static_cast<float>(itau),
                                 static_cast<float>(ix),
                                 static_cast<float>(iy),
                                 static_cast<float>(ieta),
                                 static_cast<float>(e_local*hbarc),
                                 static_cast<float>(p_local*hbarc),
                                 static_cast<float>(T_local*hbarc),
                                 static_cast<float>(ux),
                                 static_cast<float>(uy),
                                 static_cast<float>(ueta)};
                fwrite(ideal, sizeof(float), 10, out_file_xyeta);
            }
        }
    }
    fclose(out_file_xyeta);
}


