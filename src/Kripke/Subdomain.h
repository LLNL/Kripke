#ifndef KRIPKE_SUBDOMAIN_H__
#define KRIPKE_SUBDOMAIN_H__

#include <vector>

// Foreward Decl
struct Directions;
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

#endif
