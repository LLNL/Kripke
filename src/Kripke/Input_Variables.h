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
  int npx, npy, npz;            // The number of processors in x,y,z
  int nx, ny, nz;               // Number of spatial zones in x,y,z
  int num_dirsets_per_octant;
  int num_dirs_per_dirset;
  int num_groupsets;
  int num_groups_per_groupset;  //
  int num_zonesets_dim[3];      // number of zoneset in x, y, z
  int niter;                    // number of solver iterations to run
  int legendre_order;           // Scattering order (number Legendre coeff's - 1)

  Nesting_Order nesting;        // Data layout and loop ordering (of Psi)
};

#endif
