/*--------------------------------------------------------------------------
 * Header file for the Grid_Data data structures
 *--------------------------------------------------------------------------*/

#ifndef KRIPKE_GRID_DATA_H__
#define KRIPKE_GRID_DATA_H__

#include <Kripke/directions.h>
#include <Kripke/SubTVec.h>
#include <mpi.h>
#include <vector>

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

typedef std::vector< std::vector<double> > Plane_Data;

struct Input_Variables;

struct Group_Dir_Set {
  int num_groups;
  int num_directions;

  int group0;
  int direction0;

  Directions *directions;

  // Variables
  SubTVec *psi;
  SubTVec *rhs;
  SubTVec *sigt;
};

struct Grid_Data {
public:
  Grid_Data(Input_Variables *input_vars, int num_dirs, int num_grps, MPI_Comm comm);

  int num_zones;                    // Total Number of zones in this grid
  int nzones[3];                    // Number of zones in each dimension

  int mynbr[3][2];                  // Neighboring MPI ranks in each dimension

  std::vector<double> deltas[3];    // Spatial grid deltas in each dimension
  std::vector<double> volume;       // Spatial zone volumes
  
  int num_directions;
  int num_groups;

  // Group/Angle sets
  std::vector<Group_Dir_Set> gd_sets;

  // Variables:
  std::vector<double>  tmp_sigma_tot;
  
private:
  void computeGrid(int dim, int npx, int nx_g, int isub_ref, double xmin, double xmax);
};

#endif
