#ifndef KRIPKE_GRID_DATA_H__
#define KRIPKE_GRID_DATA_H__

#include <Kripke/Directions.h>
#include <Kripke/Kernel.h>
#include <Kripke/Timing.h>
#include <mpi.h>
#include <vector>

// Foreward Decl
struct Input_Variables;
struct Grid_Data;
struct SubTVec;


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

struct Neighbor{
  int mpi_rank;
  int subdomain_id;
};

/**
 * Contains parameters and variables that describe a single Group Set and
 * Direction Set.
 */
struct Subdomain {
  Subdomain();
  ~Subdomain();

  void randomizeData(void);
  void copy(Subdomain const &b);
  bool compare(Subdomain const &b, double tol, bool verbose);
  void computeSweepIndexSet(void);

  int idx_group_set;
  int idx_dir_set;
  int idx_zone_set;

  int num_groups;       // Number of groups in this set
  int num_directions;   // Number of directions in this set
  int num_zones;        // Number of zones in this set

  int nzones[3];                    // Number of zones in each dimension
  std::vector<double> deltas[3];    // Spatial grid deltas in each dimension (including ghost zones)

  int group0;           // Starting global group id
  int direction0;       // Starting global direction id

  Directions *directions;
  Grid_Sweep_Block sweep_block;

  // Neighbors
  Neighbor upwind[3];
  Neighbor downwind[3];

  // Sweep boundary data
  std::vector<double> plane_data[3];

  // Variables
  SubTVec *psi;         // Solution
  SubTVec *rhs;         // RHS, source term
  SubTVec *sigt;        // Zonal per-group cross-section

};







/**
 * Contains all grid parameters and variables.
 */
struct Grid_Data {
public:
  explicit Grid_Data(Input_Variables *input_vars);
  ~Grid_Data();

  void randomizeData(void);
  void copy(Grid_Data const &b);
  bool compare(Grid_Data const &b, double tol, bool verbose);

  Timing timing;

  int niter;

  double source_value;

  std::vector<double> sigma_tot;            // Cross section data

  int num_group_sets;                       // Number of group-sets
  int num_groups_per_set;                   // How many groups in each set
  int num_direction_sets;                   // Number of direction-sets
  int num_directions_per_set;               // Number of directions per dir set
  int num_zone_sets;

  std::vector<Directions> directions;       // Direction data

  Kernel *kernel;

  int mynbr[3][2];                  // Neighboring MPI ranks in each dimension


  // Group/Angle/Zone sets
  std::vector<Subdomain> subdomains;

  // Variables:
  int num_legendre;
  int total_num_moments;

  SubTVec *phi;               // Moments of psi
  SubTVec *phi_out;           // Scattering source (moments)

  SubTVec *ell;               // L matrix in nm_offset coordinates
  SubTVec *ell_plus;          // L+ matrix in nm_offset coordinates

private:
  std::vector<double> computeGrid(int dim, int npx, int nx_g, int isub_ref, double xmin, double xmax);
};

#endif
