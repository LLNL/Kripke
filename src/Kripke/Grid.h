#ifndef KRIPKE_GRID_DATA_H__
#define KRIPKE_GRID_DATA_H__

#include <Kripke/Directions.h>
#include <Kripke/Kernel.h>
#include <Kripke/Subdomain.h>
#include <Kripke/Timing.h>
#include <mpi.h>
#include <vector>

// Foreward Decl
struct Input_Variables;
struct Grid_Data;
struct SubTVec;


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

  // Group/Angle/Zone sets
  std::vector<Subdomain> subdomains;

  // Variables:
  int num_legendre;
  int total_num_moments;

  // Per directionset ell and ell_plus matrices (Subdomain point into these arrays)
  std::vector<SubTVec *> ell;               // L matrix in nm_offset coordinates
  std::vector<SubTVec *> ell_plus;          // L+ matrix in nm_offset coordinates

  // Per zoneset phi and phi_out (Subdomains point into these arrays)
  std::vector<SubTVec *> phi;               // Moments of psi
  std::vector<SubTVec *> phi_out;           // Scattering source
};

#endif
