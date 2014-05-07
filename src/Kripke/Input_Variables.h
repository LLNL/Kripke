/*--------------------------------------------------------------------------
 * Header file for the Input_Variables structure
 *--------------------------------------------------------------------------*/

#ifndef KRIPKE_INPUT_VARIABLES_H__
#define KRIPKE_INPUT_VARIABLES_H__

#include <Kripke.h>
#include <string>

/*--------------------------------------------------------------------------
 * Define the Input_Variables structure.
 *
 * npx                   : The number of processors in the x-direction.
 * npy                   : The number of processors in the y-direction.
 * npz                   : The number of processors in the z-direction.
 * nx                    : Number of spatial zones in the x-direction.
 * ny                    : Number of spatial zones in the y-direction.
 * nz                    : Number of spatial zones in the z-direction.
 * niter                 : The number of sweep iterations.
 * num_direction_per_octant : The number of directions per octant
 *--------------------------------------------------------------------------*/

struct Input_Variables {
  void read(std::string const &fname);
  void print(void) const;

  int npx, npy, npz;
  int nx, ny, nz;
  int num_dirsets_per_octant;
  int num_dirs_per_dirset;
  int num_groupsets;
  int num_groups_per_groupset;
  int niter;

  Nesting_Order nesting;
};

#endif
