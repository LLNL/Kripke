#ifndef KRIPKE_GRID_DATA_H__
#define KRIPKE_GRID_DATA_H__

#include <Kripke/Directions.h>
#include <Kripke/Kernel.h>
#include <mpi.h>
#include <vector>

// Foreward Decl
struct Input_Variables;
struct Grid_Data;
struct SubTVec;


/**
 * Contains parameters and variables that describe a single Group Set and
 * Direction Set.
 */
struct Group_Dir_Set {
  Group_Dir_Set();
  ~Group_Dir_Set();

  void allocate(Grid_Data *grid_data, Nesting_Order nesting);
  void randomizeData(void);
  void copy(Group_Dir_Set const &b);
  bool compare(int gs, int ds, Group_Dir_Set const &b, double tol, bool verbose);

  int num_groups;       // Number of groups in this set
  int num_directions;   // Number of directions in this set

  int group0;           // Starting global group id
  int direction0;       // Starting global direction id

  Directions *directions;

  // Variables
  SubTVec *psi;         // Solution
  SubTVec *rhs;         // RHS, source term
};



/**
 * Provides sweep index sets for a given octant.
 * This generalizes the sweep pattern, and allows for experimenting with
 * a tiled approach to on-node sweeps.
 */
struct Grid_Sweep_Block {
  int start_i, start_j, start_k; // starting index
  int end_i, end_j, end_k; // termination conditon (one past)
  int inc_i, inc_j, inc_k; // increment
};



/**
 * Contains all grid parameters and variables.
 */
struct Grid_Data {
public:
  Grid_Data(Input_Variables *input_vars, Directions *directions);
  ~Grid_Data();

  void randomizeData(void);
  void copy(Grid_Data const &b);
  bool compare(Grid_Data const &b, double tol, bool verbose);

  int num_zones;                    // Total Number of zones in this grid
  int nzones[3];                    // Number of zones in each dimension

  int mynbr[3][2];                  // Neighboring MPI ranks in each dimension

  std::vector<double> deltas[3];    // Spatial grid deltas in each dimension

  // Sweep index sets for each octant
  std::vector<Grid_Sweep_Block> octant_extent;

  // Group/Angle sets
  std::vector< std::vector<Group_Dir_Set> > gd_sets;

  // Variables:
  int num_moments;
  SubTVec *sigt;              // Zonal per-group cross-section
  SubTVec *phi;               // Moments of psi
  SubTVec *phi_out;           // Scattering source (moments)
  SubTVec *ell;               // L matrix in nm_offset coordinates
  SubTVec *ell_plus;          // L+ matrix in nm_offset coordinates
  std::vector<int> nm_table;  // n, m indicies for traversing ell, ell_plus

private:
  void computeGrid(int dim, int npx, int nx_g, int isub_ref, double xmin, double xmax);
  void computeSweepIndexSets(int block_size);
};

#endif
