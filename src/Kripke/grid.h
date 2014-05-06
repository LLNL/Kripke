/*--------------------------------------------------------------------------
 * Header file for the Grid_Data data structures
 *--------------------------------------------------------------------------*/

#ifndef KRIPKE_GRID_DATA_H__
#define KRIPKE_GRID_DATA_H__

#include <Kripke/directions.h>
#include <Kripke/Kernel.h>
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

struct Input_Variables;
struct Grid_Data;

struct SubTVec;
struct LMat;

struct Group_Dir_Set {
  Group_Dir_Set();
  ~Group_Dir_Set();

  void allocate(Grid_Data *grid_data, Nesting_Order nesting);

  int num_groups;
  int num_directions;

  int group0;
  int direction0;

  Directions *directions;

  // Variables
  SubTVec *psi;         // Solution
  SubTVec *rhs;         // RHS, source term
  SubTVec *sigt;        // Zonal per-group cross-section
};

struct Grid_Data {
public:
  Grid_Data(Input_Variables *input_vars, Directions *directions);
  ~Grid_Data();

  int num_zones;                    // Total Number of zones in this grid
  int nzones[3];                    // Number of zones in each dimension

  int mynbr[3][2];                  // Neighboring MPI ranks in each dimension

  std::vector<double> deltas[3];    // Spatial grid deltas in each dimension
  std::vector<double> volume;       // Spatial zone volumes

  // Group/Angle sets
  std::vector< std::vector<Group_Dir_Set> > gd_sets;

  // Variables:
  int num_moments;
  SubTVec *phi;         // Moments of psi
  LMat *ell;         // L matrix
  LMat *ell_plus;    // L+ matrix
  
private:
  void computeGrid(int dim, int npx, int nx_g, int isub_ref, double xmin, double xmax);
};

#endif
