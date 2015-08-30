/*--------------------------------------------------------------------------
 * Header file for the Input_Variables structure
 *--------------------------------------------------------------------------*/

#ifndef KRIPKE_INPUT_VARIABLES_H__
#define KRIPKE_INPUT_VARIABLES_H__

#include<Kripke.h>

/**
 * This structure defines the input parameters to setup a problem.
 */

struct Input_Variables {
  Input_Variables();
  
  bool checkValues(void) const;
  
  // Problem Description
  int nx, ny, nz;               // Number of spatial zones in x,y,z
  int num_directions;           // Total number of directions
  int num_groups;               // Total number of energy groups
  int legendre_order;           // Scattering order (number Legendre coeff's - 1)
  int quad_num_polar;           // Number of polar quadrature points (0 for dummy)
  int quad_num_azimuthal;       // Number of azimuthal quadrature points (0 for dummy)

  // On-Node Options
  Nesting_Order nesting;        // Data layout and loop ordering (of Psi)
  
  // Parallel Decomp
  int npx, npy, npz;            // The number of processors in x,y,z
  int num_dirsets;              // Number of direction sets
  int num_groupsets;            // Number of energy group sets
  int num_zonesets_dim[3];      // Number of zoneset in x, y, z  
  int layout_pattern;           // Which subdomain/task layout to use
  
  // Physics and Solver Options
  int niter;                    // number of solver iterations to run
  ParallelMethod parallel_method;
  double sigt[3];               // total cross section for 3 materials
  double sigs[3];               // total scattering cross section for 3 materials
  
  // Output Options
  std::string run_name;         // Name to use when generating output files
  std::string outfile;          // name of output file
#ifdef KRIPKE_USE_SILO
  std::string silo_basename;    // name prefix for silo output files
#endif

};

#endif
