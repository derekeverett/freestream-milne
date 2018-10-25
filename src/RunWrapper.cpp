#include "FreestreamMilne.cpp"

int main(void)
{
  //Declare an instance of FREESTREAMMILNE class
  FREESTREAMMILNE fsmilne;

  //try to initialize the energy density from a vector of DIMX * DIMY * DIMETA elements
  /*
  int npoints = 101 * 101;
  size_t size = npoints * sizeof(double);
  std::vector<double> init_e(size); //initial energy density vector
  std::vector<double> final_e(size); //final energy density vector
  std::vector<double> final_p(size); //final pressure vector
  std::vector<double> final_ut(size); // etc...
  std::vector<double> final_ux(size);
  std::vector<double> final_uy(size);
  std::vector<double> final_un(size);
  std::vector<double> final_pitt(size);
  std::vector<double> final_pitx(size);
  std::vector<double> final_pity(size);
  std::vector<double> final_pitn(size);
  std::vector<double> final_pixx(size);
  std::vector<double> final_pixy(size);
  std::vector<double> final_pixn(size);
  std::vector<double> final_piyy(size);
  std::vector<double> final_piyn(size);
  std::vector<double> final_pinn(size);
  std::vector<double> final_Pi(size);

  //this step would be handled by another module... e.g. TRENTO
  for (int i = 0; i < npoints; i++) init_e[i] = 1.0; //initialize to nontrivial values


  //tell fsmilne how big the grid is
  fsmilne.gridSize = npoints;

  //pass the initial energy density vector
  //this function is only used when IC.ENERGY == 5 ,the option for an initial energy density from vector
  fsmilne.initialize_from_vector(init_e);
  */
  //run the freestreaming evolution
  fsmilne.run_freestream_milne();

  /*
  //grab the final hydro vectors to pass to another module
  fsmilne.output_to_vectors(final_e,
                            final_p,
                            final_ut,
                            final_ux,
                            final_uy,
                            final_un,
                            final_pitt,
                            final_pitx,
                            final_pity,
                            final_pitn,
                            final_pixx,
                            final_pixy,
                            final_pixn,
                            final_piyy,
                            final_piyn,
                            final_pinn,
                            final_Pi);
  */
  printf("run_freestream_milne() ran sucessfully \n");

  /*
  //check that the vectors were upated
  for (int i = 0; i < npoints; i++)
  {
    printf("e [ %d ] = %f \n", i, final_e[i]);
  }
  */
}
