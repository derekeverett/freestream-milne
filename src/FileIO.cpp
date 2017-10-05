#include <unistd.h>
#include <stdio.h>
#include <fstream>

void writeScalarToFile(float *var, char name[255])
{
  char cwd[1024];
  getcwd(cwd, sizeof(cwd));
  std::ofstream myfile;
  char filename[255] = "";
  sprintf(filename, "output/%s.dat", name);
  myfile.open(filename);
  for (int ieta = 0; ieta < DIM_ETA; ieta++)
  {
    for (int iy = 0; iy < DIM_Y; iy++)
    {
      for (int ix = 0; ix < DIM_X; ix++)
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

void writeVectorToFile(float **var, char name[255], int idx)
{
  char cwd[1024];
  getcwd(cwd, sizeof(cwd));
  std::ofstream myfile;
  char filename[255] = "";
  sprintf(filename, "output/%s.dat", name);
  myfile.open(filename);
  for (int ieta = 0; ieta < DIM_ETA; ieta++)
  {
    for (int iy = 0; iy < DIM_Y; iy++)
    {
      for (int ix = 0; ix < DIM_X; ix++)
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
void writeScalarToFileProjection(float *var, char name[255])
{
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
      int ieta = (DIM_ETA - 1) / 2; // at eta = 0
      int is = (DIM_Y * DIM_ETA) * ix + (DIM_ETA) * iy + ieta; //the column packed index spanning x, y, eta
      myfile << var[is] << " "; //different columns for y values
    }
    myfile << "\n"; // different rows correspond to different x values
  }
  myfile.close();
}

void writeVectorToFileProjection(float **var, char name[255], int idx)
{
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
      int ieta = (DIM_ETA - 1) / 2; //at eta = 0
      int is = (DIM_Y * DIM_ETA) * ix + (DIM_ETA) * iy + ieta; //the column packed index spanning x, y, z
      myfile << var[idx][is] << " "; //different columns for y values
    }
    myfile << "\n"; // different rows correspond to different x values
  }
  myfile.close();
}
