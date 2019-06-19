freestream-milne (c) Derek Everett

## Purpose
freestream-milne numerically solves the collisionless Boltzmann Equation assuming that
1. the particles are massless
2. initial distribution function is isotropic in azimuthal momentum space $px, py$
3. the initial distribution is a gaussian $(y - \eta)$ in 3D case, or delta function in $(y-\eta)$ in 2D case. 

Given an initial energy density profile $T^{00}$, it computes the stress energy tensor $T^{\mu\nu}$ at a later time, and also a viscous hydrodynamic decomposition. 

OpenMP is used to accelerate the computation.

## Usage 

Parameters and options can be set in freestream_input, which is read in at run time.
They are described below
OUTPUTFORMAT (int) 1 : disk , 2 : write to vectors in memory
BARYON (int) 0 : freestream the energy density, 1 : freestream the baryon density and energy density
IC_ENERGY (int) 1 : Smooth Gaussian, 2 : MC Guassian, 3: read tab delimited file for energy density, 4 : read superMC block format, 5 : read from vector in Memory
IC_BARYON (int) 1 : Smooth Gaussian, 2 : MC Guassian, 3: read tab delimited file for baryon density
ETA_WIDTH (float) used to create a 3D profile of energy density from 2D MC Glb (only for superMC block)
ETA_FLAT (float) used to create a 3D profile of energy density from 2D MC Glb (only for superMC block)
SIGMA (float) std. deviation of one particle distribution function f ~ exp((y - eta)^2 / 2 SIGMA^2)
SIGMA_B (float) std. deviation of one particle distribution function of baryon charges f_B ~ exp((y - eta)^2 / 2 SIGMA_B^2)
DIM_X (int) number of grid points in x
DIM_Y (int) number of grid points in y
DIM_ETA (int) number of grid points in eta
DIM_RAP (int) number of grid points in rapidity for riemann sum
DIM_PHIP (int) number of grid points in phi_p for riemann sum
DX (float) lattice spacing in x [fm]
DY (float) lattice spacing in y [fm]
DETA (float) lattice spacing in eta
DRAP (float) lattice spacing in rapidity (differential) DEPRECATED
DTAU (float) free streaming proper time interval [fm/c]
TAU0 (float) start of free streaming proper time [fm/c]
EOS_TYPE (int) 1 : conformal EoS, 2 : Wuppertal-Bhudapest Parameterization
E_FREEZE (float) [GeV/fm^3] used to calculate fraction of energy contained within surface of energy density > E_FREEZE
