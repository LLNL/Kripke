/*--------------------------------------------------------------------------
 * Header file for the Grid_Data data structures
 *--------------------------------------------------------------------------*/

#ifndef included_grid_data
#define included_grid_data

#include "directions.h"

/*--------------------------------------------------------------------------
 * Define the Grid_Data structure.
 *
 * nzones         : The number of zones used in each coordinate direction
 * nprocs         : The number of processors used in each coordinate direction
 * mynbr          : MPI process id's of neighbors for each grid face
 * coord_lo       : Local x,y,z lower coordinates
 * coord_hi       : Local x,y,z upper coordinates
 * xmin           : The miniumum spatial value for the x direction
 * xmax           : The maximum spatial value for the x direction
 * ymin           : The miniumum spatial value for the y direction
 * ymax           : The maximum spatial value for the y direction
 * zmin           : The miniumum spatial value for the z direction
 * zmax           : The maximum spatial value for the z direction
 * deltas         : The delta_x, delta_y, and delta_z arrays
 * volume         : The volume of each spatial zone
 * eps            : Unit Roundoff
 * global_num_zones: The global number of spatial zones
 * num_directions : The number of directions
 * ilower         : Lower grid indices
 *--------------------------------------------------------------------------*/

typedef struct {
  int    *nzones;
  int    *nprocs;

  int mynbr[3][2];

  double coord_lo[3];
  double coord_hi[3];
  double xmin, xmax, ymin, ymax, zmin, zmax;

  double *deltas[3];
  double *volume;
  double eps;   /* unit roundoff */

  double global_num_zones;
  int num_directions;
  int ilower[3];
  int iupper[3];

  Directions *directions;
}   Grid_Data;

#endif
