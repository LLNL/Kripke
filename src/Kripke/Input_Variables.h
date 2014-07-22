/*--------------------------------------------------------------------------
 * Header file for the Input_Variables structure
 *--------------------------------------------------------------------------*/

#ifndef KRIPKE_INPUT_VARIABLES_H__
#define KRIPKE_INPUT_VARIABLES_H__

#include <Kripke.h>
#include <string>

/**
 * This structure defines the input parameters to setup a problem.
 *
 * npx                   :
 * npy                   : The number of processors in the y-direction.
 * npz                   : The number of processors in the z-direction.
 * nx                    : Number of spatial zones in the x-direction.
 * ny                    : Number of spatial zones in the y-direction.
 * nz                    : Number of spatial zones in the z-direction.
 * niter                 : The number of sweep iterations.
 * num_direction_per_octant : The number of directions per octant
 */

struct Input_Variables {
  int npx, npy, npz;            // The number of processors in x,y,z
  int nx, ny, nz;               // Number of spatial zones in x,y,z
  int num_dirsets_per_octant;
  int num_dirs_per_dirset;
  int num_groupsets;
  int num_groups_per_groupset;
  int niter;                    // number of solver iterations to run
  int legendre_order;

  Nesting_Order nesting;
  int block_size;               // Experimental: used for tiling spatial zones
};

#endif
