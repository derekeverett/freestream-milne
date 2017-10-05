//trigTable is a table with 10 rows for ten combinations or p_(mu)p_(nu) normalized by the energy
#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_sort_vector.h>

#define REGULATE 1 // 1 to regulate dilute regions of space, sets energy density to zero if < tolerance
void calculateTrigTable(float ***trigTable)
{
  for (int ithetap = 0; ithetap < DIM_THETAP; ithetap++)
  {
    for (int iphip = 0; iphip < DIM_PHIP; iphip++)
    {
      float thetap = float(ithetap) * (PI) / float(DIM_THETAP);
      float phip = float(iphip) * (2.0 * PI) / float(DIM_PHIP);
      trigTable[0][ithetap][iphip] = 1.0; //p_t, p_t component
      trigTable[1][ithetap][iphip] = sin(thetap) * cos(phip); //p_t, p_x
      trigTable[2][ithetap][iphip] = sin(thetap) * sin(phip); //p_t, p_y
      trigTable[3][ithetap][iphip] = cos(thetap); //p_t, p_z
      trigTable[4][ithetap][iphip] = sin(thetap) * cos(phip) * sin(thetap) * cos(phip); //p_x, p_x
      trigTable[5][ithetap][iphip] = sin(thetap) * cos(phip) * sin(thetap) * sin(phip); //p_x, p_y
      trigTable[6][ithetap][iphip] = sin(thetap) * cos(phip) * cos(thetap); //p_x, p_z
      trigTable[7][ithetap][iphip] = sin(thetap) * sin(phip) * sin(thetap) * sin(phip); //p_y, p_y
      trigTable[8][ithetap][iphip] = sin(thetap) * sin(phip) * cos(thetap); //p_y, p_z
      trigTable[9][ithetap][iphip] = cos(thetap) * cos(thetap); //p_z, p_z
      trigTable[10][ithetap][iphip] = sin(thetap); //just sin (thetap), useful for when we calculate stress tensor
    }
  }
}
void calculateStressTensor(float **stressTensor, float ***shiftedDensity, float ***trigTable)
{
  float d_phip = (2.0 * PI) / float(DIM_PHIP);
  for (int ivar = 0; ivar < 10; ivar++) //ten independent components
  {
    #pragma omp parallel for simd
    for (int is = 0; is < DIM; is++) //the column packed index for x, y and z
    {
      for (int irap = 0; irap < DIM_RAP; irap++)
      {
        for (int iphip = 0; iphip < DIM_PHIP; iphip++)
        {
          //rather than gauss quadrature, just doing a elementary Riemann sum here; check convergence!
          stressTensor[ivar][is] += shiftedDensity[is][irap][iphip] * hypertrigTable[ivar][irap][iphip] * hypertrigTable[10][irap][iphip] * DRAP * d_phip;
        }
      }
    }
  }
}
void solveEigenSystem(float **stressTensor, float *energyDensity, float **flowVelocity)
{
  float tolerance = 1e18; //set quantities to zero which are less than 10^(-18)

  #pragma omp parallel for simd
  for (int is = 0; is < DIM; is++)
  {
    gsl_matrix *Tmunu; //T^(mu,nu) with two contravariant indices; we need to lower an index
    //using the metric to find the eigenvectors of T^(mu)_(nu) with one contravariant and one contravariant index
    Tmunu = gsl_matrix_alloc(4,4);
    gsl_matrix *gmunu;
    gmunu = gsl_matrix_alloc(4,4);
    gsl_matrix_complex *eigen_vectors;
    eigen_vectors = gsl_matrix_complex_alloc(4,4);
    gsl_vector_complex *eigen_values;
    eigen_values = gsl_vector_complex_alloc(4);
    //set the values of the energy momentum tensor
    gsl_matrix_set(Tmunu, 0, 0, stressTensor[0][is]); //tt
    gsl_matrix_set(Tmunu, 0, 1, stressTensor[1][is]); //tx
    gsl_matrix_set(Tmunu, 0, 2, stressTensor[2][is]); //ty
    gsl_matrix_set(Tmunu, 0, 3, stressTensor[3][is]); //tz
    gsl_matrix_set(Tmunu, 1, 1, stressTensor[4][is]); //xx
    gsl_matrix_set(Tmunu, 1, 2, stressTensor[5][is]); //xy
    gsl_matrix_set(Tmunu, 1, 3, stressTensor[6][is]); //xz
    gsl_matrix_set(Tmunu, 2, 2, stressTensor[7][is]); //yy
    gsl_matrix_set(Tmunu, 2, 3, stressTensor[8][is]); //yz
    gsl_matrix_set(Tmunu, 3, 3, stressTensor[9][is]); //zz
    gsl_matrix_set(Tmunu, 1, 0, stressTensor[1][is]); //xt
    gsl_matrix_set(Tmunu, 2, 0, stressTensor[2][is]); //yt
    gsl_matrix_set(Tmunu, 3, 0, stressTensor[3][is]); //zt
    gsl_matrix_set(Tmunu, 2, 1, stressTensor[5][is]); //yx
    gsl_matrix_set(Tmunu, 3, 1, stressTensor[6][is]); //zx
    gsl_matrix_set(Tmunu, 3, 2, stressTensor[8][is]); //zy

    //set the values of the "metric"; not really the metric, but the numerical constants
    //which are multiplied by the elements of T^(mu,nu) to get the values of T^(mu)_(nu)
    gsl_matrix_set(gmunu, 0, 0, 1.0); //tt
    gsl_matrix_set(gmunu, 0, 1, -1.0); //tx
    gsl_matrix_set(gmunu, 0, 2, -1.0); //ty
    gsl_matrix_set(gmunu, 0, 3, -1.0); //tz
    gsl_matrix_set(gmunu, 1, 0, 1.0); //xt
    gsl_matrix_set(gmunu, 1, 1, -1.0); //xx
    gsl_matrix_set(gmunu, 1, 2, -1.0); //xy
    gsl_matrix_set(gmunu, 1, 3, -1.0); //xz
    gsl_matrix_set(gmunu, 2, 0, 1.0); //yt
    gsl_matrix_set(gmunu, 2, 1, -1.0); //yx
    gsl_matrix_set(gmunu, 2, 2, -1.0); //yy
    gsl_matrix_set(gmunu, 2, 3, -1.0); //yz
    gsl_matrix_set(gmunu, 3, 0, 1.0); //zt
    gsl_matrix_set(gmunu, 3, 1, -1.0); //zx
    gsl_matrix_set(gmunu, 3, 2, -1.0); //zy
    gsl_matrix_set(gmunu, 3, 3, -1.0); //zz
    //lower one index of the stress tensor; save it to the same matrix to save memory
    gsl_matrix_mul_elements(Tmunu, gmunu); //result stored in Tmunu !this multiplies element-wise, not ordinary matrix multiplication!
    gsl_eigen_nonsymmv_workspace *eigen_workspace;
    eigen_workspace = gsl_eigen_nonsymmv_alloc(4);
    gsl_eigen_nonsymmv(Tmunu, eigen_values, eigen_vectors, eigen_workspace);
    gsl_eigen_nonsymmv_free(eigen_workspace);

    //***does this have a solution for energy density and flow at every point?
    for (int i = 0; i < 4; i++)
    {
      gsl_complex eigenvalue = gsl_vector_complex_get(eigen_values, i);

      if (GSL_REAL(eigenvalue) > 0.0 && GSL_IMAG(eigenvalue) == 0.0)
      {
        double v0 = GSL_REAL(gsl_matrix_complex_get(eigen_vectors, i , 0));
        double v1 = GSL_REAL(gsl_matrix_complex_get(eigen_vectors, i , 1));
        double v2 = GSL_REAL(gsl_matrix_complex_get(eigen_vectors, i , 2));
        double v3 = GSL_REAL(gsl_matrix_complex_get(eigen_vectors, i , 3));
        //double euclideanLength = v0*v0 + v1*v1 + v2*v2 + v3*v3; //gsl normalizes eigenvectors to euclideanLength = 1; this is just a check
        double minkowskiLength = v0*v0 - (v1*v1 + v2*v2 + v3*v3); //we want to flow velocity normalized s.t. minkowskiLength = 1
        double scaleFactor = 1.0 / sqrt(minkowskiLength); //so we need to scale all the elements of the eigenvector by scaleFactor
        //printf("scaled eigenvector %d is (%f ,%f , %f, %f) and eigenvalue %d is %f\n", i, v0, v1, v2, v3, i, GSL_REAL(eigenvalue));
        //set values of energy density and flow velocity

        if ((REGULATE) && (GSL_REAL(eigenvalue) * tolerance <= 1.0)) //regulate dilute regions
        {
          energyDensity[is] = 0.0;
          flowVelocity[0][is] = 0.0;
          flowVelocity[1][is] = 0.0;
          flowVelocity[2][is] = 0.0;
          flowVelocity[3][is] = 0.0;
        }
        energyDensity[is] = GSL_REAL(eigenvalue) / scaleFactor; //do we need to scale the eigenvalue by the inverse of scaleFactor?
        flowVelocity[0][is] = v0 * scaleFactor;
        flowVelocity[1][is] = v1 * scaleFactor;
        flowVelocity[2][is] = v2 * scaleFactor;
        flowVelocity[3][is] = v3 * scaleFactor;
      }
    }
  }
}
void calculateBulkPressure(float **stressTensor, float *energyDensity, float *pressure, float *bulkPressure)
{
  for (int is = 0; is < DIM; is++)
  {
    // PI = -1/3 * (T^(mu)_(mu) - epsilon) - p
    // T^(mu)_(mu) = T^(0,0) - T^(1,1) - T^(2,2) - T^(3,3)
    float a =  stressTensor[0][is] - stressTensor[4][is] - stressTensor[7][is] - stressTensor[9][is];
    bulkPressure[is] = (-1.0/3.0) * (a - energyDensity[is]) - pressure[is];
  }
}
void calculateShearViscTensor(float **stressTensor, float *energyDensity, float **flowVelocity, float *pressure, float *bulkPressure, float **shearTensor)
{
  for (int is = 0; is < DIM; is++)
  {
    // pi^(mu,nu) = T^(mu,nu) - epsilon * u^(mu)u^(nu) + (P + PI) * (g^(mu,nu) - u^(mu)u^(nu))
    float c = energyDensity[is] + pressure[is] + bulkPressure[is];
    shearTensor[0][is] = stressTensor[1][is] - flowVelocity[0][is] * flowVelocity[1][is] * c; //pi^(0,1)
    shearTensor[1][is] = stressTensor[2][is] - flowVelocity[0][is] * flowVelocity[2][is] * c; //pi^(0,2)
    shearTensor[2][is] = stressTensor[3][is] - flowVelocity[0][is] * flowVelocity[3][is] * c; //pi^(0,3)
    shearTensor[3][is] = stressTensor[5][is] - flowVelocity[1][is] * flowVelocity[2][is] * c; //pi^(1,2)
    shearTensor[4][is] = stressTensor[6][is] - flowVelocity[1][is] * flowVelocity[3][is] * c; //pi^(1,3)
    shearTensor[5][is] = stressTensor[8][is] - flowVelocity[2][is] * flowVelocity[3][is] * c; //pi^(2,3)
  }
}
