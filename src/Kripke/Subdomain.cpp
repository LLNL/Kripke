#include <Kripke/Grid.h>
#include <Kripke/SubTVec.h>
#include <Kripke/Input_Variables.h>

#include <cmath>
#include <sstream>

Subdomain::Subdomain() :
  idx_dir_set(0),
  idx_group_set(0),
  idx_zone_set(0),
  num_groups(0),
  num_directions(0),
  group0(0),
  direction0(0),
  directions(NULL),
  psi(NULL),
  rhs(NULL),
  sigt(NULL)
{
  plane_data[0] = NULL;
  plane_data[1] = NULL;
  plane_data[2] = NULL;
}
Subdomain::~Subdomain(){
  delete psi;
  delete rhs;
  delete sigt;
  delete plane_data[0];
  delete plane_data[1];
  delete plane_data[2];
}


/**
  Setup subdomain and allocate data
*/
void Subdomain::setup(int sdom_id, Input_Variables *input_vars, int gs, int ds, int zs, std::vector<Directions> &direction_list, Kernel *kernel, int mynbr[3][2]){
  // set the set indices
  idx_group_set = gs;
  idx_dir_set = ds;
  idx_zone_set = zs;

  num_groups = input_vars->num_groups_per_groupset;
  group0 = gs * input_vars->num_groups_per_groupset;

  num_directions = input_vars->num_dirs_per_dirset;
  direction0 = ds * input_vars->num_dirs_per_dirset;
  directions = &direction_list[direction0];

  nzones[0] = deltas[0].size()-2;
  nzones[1] = deltas[1].size()-2;
  nzones[2] = deltas[2].size()-2;
  num_zones = nzones[0] * nzones[1] * nzones[2];

  // allocate storage for the sweep boundary data
  plane_data[0] = new SubTVec(kernel->nestingPsi(), num_groups, num_directions, nzones[1] * nzones[2]);
  plane_data[1] = new SubTVec(kernel->nestingPsi(), num_groups, num_directions, nzones[0] * nzones[2]);
  plane_data[2] = new SubTVec(kernel->nestingPsi(), num_groups, num_directions, nzones[0] * nzones[1]);

  // allocate the storage for solution and source terms
  psi = new SubTVec(kernel->nestingPsi(), num_groups, num_directions, num_zones);
  rhs = new SubTVec(kernel->nestingPsi(), num_groups, num_directions, num_zones);
  sigt = new SubTVec(kernel->nestingSigt(), num_groups, 1, num_zones);

  computeSweepIndexSet();

  // Setup neighbor data
  int dirs[3] = { directions[0].id, directions[0].jd, directions[0].kd};
  for(int dim = 0;dim < 3;++ dim){
    if(dirs[dim] > 0){
      downwind[dim].mpi_rank = mynbr[dim][1];
      downwind[dim].subdomain_id = sdom_id;
      upwind[dim].mpi_rank = mynbr[dim][0];
      upwind[dim].subdomain_id = sdom_id;
    }
    else{
      downwind[dim].mpi_rank = mynbr[dim][0];
      downwind[dim].subdomain_id = sdom_id;
      upwind[dim].mpi_rank = mynbr[dim][1];
      upwind[dim].subdomain_id = sdom_id;
    }
  }
}



/**
 * Randomizes data for a set.
 */
void Subdomain::randomizeData(void){
  psi->randomizeData();
  rhs->randomizeData();
  sigt->randomizeData();

  for(int d = 0;d < 3;++ d){
    for(int i = 0;i < deltas[d].size();++ i){
      deltas[d][i] = drand48();
    }
  }
}

/**
 * Copies two sets, allowing for different nestings.
 */
void Subdomain::copy(Subdomain const &b){
  psi->copy(*b.psi);
  rhs->copy(*b.rhs);
  sigt->copy(*b.sigt);

  for(int d = 0;d < 3;++ d){
    deltas[d] = b.deltas[d];
  }
}

/**
 * Compares two sets, allowing for different nestings.
 */
bool Subdomain::compare(Subdomain const &b, double tol, bool verbose){
  std::stringstream namess;
  namess << "gdset[gs=" << idx_group_set << ", ds=" << idx_dir_set << ", zs=" << idx_zone_set << "]";
  std::string name = namess.str();

  bool is_diff = false;
  is_diff |= psi->compare(name+".psi", *b.psi, tol, verbose);
  is_diff |= rhs->compare(name+".rhs", *b.rhs, tol, verbose);
  is_diff |= sigt->compare(name+".sigt", *b.sigt, tol, verbose);

  is_diff |= compareVector(name+".deltas[0]", deltas[0], b.deltas[0], tol, verbose);
  is_diff |= compareVector(name+".deltas[1]", deltas[1], b.deltas[1], tol, verbose);
  is_diff |= compareVector(name+".deltas[2]", deltas[2], b.deltas[2], tol, verbose);

  return is_diff;
}

/**
 * Compute sweep index sets.
 * Determines logical indices, and increments for i,j,k based on grid
 * information and quadrature set sweeping direction.
 */
void Subdomain::computeSweepIndexSet(void){
  if(directions[0].id > 0){
    sweep_block.start_i = 0;
    sweep_block.end_i = nzones[0]-1;
    sweep_block.inc_i = 1;
  }
  else {
    sweep_block.start_i = nzones[0]-1;
    sweep_block.end_i = 0;
    sweep_block.inc_i = -1;
  }

  if(directions[0].jd > 0){
    sweep_block.start_j = 0;
    sweep_block.end_j = nzones[1]-1;
    sweep_block.inc_j = 1;
  }
  else {
    sweep_block.start_j = nzones[1]-1;
    sweep_block.end_j = 0;
    sweep_block.inc_j = -1;
  }

  if(directions[0].kd > 0){
    sweep_block.start_k = 0;
    sweep_block.end_k = nzones[2]-1;
    sweep_block.inc_k =  1;
  }
  else {
    sweep_block.start_k = nzones[2]-1;
    sweep_block.end_k = 0;
    sweep_block.inc_k = -1;
  }
}



