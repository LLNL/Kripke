/*--------------------------------------------------------------------------
 * Header file for the Input_Variables structure
 *--------------------------------------------------------------------------*/

#ifndef KRIPKE_INPUT_VARIABLES_H__
#define KRIPKE_INPUT_VARIABLES_H__

#include<string>

/*--------------------------------------------------------------------------
 * Define the Input_Variables structure.
 *
 * npx                   : The number of processors in the x-direction.
 * npy                   : The number of processors in the y-direction.
 * npz                   : The number of processors in the z-direction.
 * xmin                  : Left domain value of x in physical space.
 * xmax                  : Right domain value of x in physical space.
 * ymin                  : Front domain value of y in physical space.
 * ymax                  : Back domain value of y in physical space.
 * zmin                  : Lower domain value of z in physical space.
 * zmax                  : Upper domain value of z in physical space.
 * nx                    : Number of spatial zones in the x-direction.
 * ny                    : Number of spatial zones in the y-direction.
 * nz                    : Number of spatial zones in the z-direction.
 * nlevels_kba           : If > 0, then forces the DD sweep in 3D to use
 *                         KBA for sweeps with nlevels_kba levels.
 * ncalls                : The number of calls to the sweep driver routine.
 * num_direction_per_octant : The number of directions per octant
 * bndry_types           : Boundary condition types for faces.
 *                            Type = 0  => Dirichlet
 *                            Type = 1  => Reflecting
 * bndry_values          : Boundary condition values for all groups and faces.
 * source_value          : The source value (assumed constant over
 *                         the domain) for the Boltzmann equation
 * sigma_total_value     : Value for sigma_total constant
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
  int ncalls;

  double xmin, xmax, ymin, ymax, zmin, zmax;
  double source_value;
  double sigma_total_value;
  double bndry_values[6];
  int bndry_types[6];

  char run_name[256];

};

#endif
