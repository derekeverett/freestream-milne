#include "FreestreamMilne.cpp"

int main(void)
{
  //Declare an instance of FREESTREAMMILNE class
  FREESTREAMMILNE free;

  //try to initialize the energy density from a vector of DIMX * DIMY * DIMETA elements
  int npoints = 101 * 101;
  size_t size = npoints * sizeof(float);
  std::vector<float> init(size);
  for (int i = 0; i < npoints; i++) init[i] = 1.0;

  //pass the initial energy density vector
  //this function is only used when IC.ENERGY == 5 ,the option for an initial energy density from vector
  free.initialize_from_vector(init);

  //run the freestreaming evolution
  free.run_freestream_milne();

}
