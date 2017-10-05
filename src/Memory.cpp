float ** calloc2dArray(float **array, int dim1, int dim2)
{
  array = (float **)calloc(dim1, sizeof(float *));
  for (int i1 = 0; i1 < dim1; i1++)
  {
    array[i1] = (float *)calloc(dim2, sizeof(float));
  }
  return array;
}

float *** calloc3dArray(float ***array, int dim1, int dim2, int dim3)
{
  array = (float ***)calloc(dim1, sizeof(float **));
  for (int i1 = 0; i1 < dim1; i1++)
  {
    array[i1] = (float **)calloc(dim2, sizeof(float *));
    for (int i2 = 0; i2 < dim2; i2++)
    {
      array[i1][i2] = (float *)calloc(dim3, sizeof(float));
    }
  }
  return array;
}

float **** calloc4dArray(float ****array, int dim1, int dim2, int dim3, int dim4)
{
  array = (float ****)calloc(dim1, sizeof(float ***));
  for (int i1 = 0; i1 < dim1; i1++)
  {
    array[i1] = (float ***)calloc(dim2, sizeof(float **));
    for (int i2 = 0; i2 < dim2; i2++)
    {
      array[i1][i2] = (float **)calloc(dim3, sizeof(float *));
      for (int i3 = 0; i3 < dim3; i3++)
      {
        array[i1][i2][i3] = (float *)calloc(dim4, sizeof(float));
      }
    }
  }
  return array;
}

float ***** calloc5dArray(float *****array, int dim1, int dim2, int dim3, int dim4, int dim5)
{
  array = (float *****)calloc(dim1, sizeof(float ****));
  for (int i1 = 0; i1 < dim1; i1++)
  {
    array[i1] = (float ****)calloc(dim2, sizeof(float ***));
    for (int i2 = 0; i2 < dim2; i2++)
    {
      array[i1][i2] = (float ***)calloc(dim3, sizeof(float **));
      for (int i3 = 0; i3 < dim3; i3++)
      {
        array[i1][i2][i3] = (float **)calloc(dim4, sizeof(float *));
        for (int i4 = 0; i4 < dim4; i4++)
        {
          array[i1][i2][i3][i4] = (float *)calloc(dim5, sizeof(float));
        }
      }
    }
  }
  return array;
}

void free2dArray(float **array, int dim1)
{
  for (int i1 = 0; i1 < dim1; i1++)
  {
    free(array[i1]);
  }
  free(array);
}

void free3dArray(float ***array, int dim1, int dim2)
{
  for (int i1 = 0; i1 < dim1; i1++)
  {
    for (int i2 = 0; i2 < dim2; i2++)
    {
      free(array[i1][i2]);
    }
    free(array[i1]);
  }
  free(array);
}
void free4dArray(float ****array, int dim1, int dim2, int dim3)
{
  for (int i1 = 0; i1 < dim1; i1++)
  {
    for (int i2 = 0; i2 < dim2; i2++)
    {
      for (int i3 = 0; i3 < dim3; i3++)
      {
        free(array[i1][i2][i3]);
      }
      free(array[i1][i2]);
    }
    free(array[i1]);
  }
  free(array);
}

void free5dArray(float *****array, int dim1, int dim2, int dim3, int dim4)
{
  for (int i1 = 0; i1 < dim1; i1++)
  {
    for (int i2 = 0; i2 < dim2; i2++)
    {
      for (int i3 = 0; i3 < dim3; i3++)
      {
        for (int i4 = 0; i4 < dim4; i4++)
        {
          free(array[i1][i2][i3][i4]);
        }
        free(array[i1][i2][i3]);
      }
      free(array[i1][i2]);
    }
    free(array[i1]);
  }
  free(array);
}
