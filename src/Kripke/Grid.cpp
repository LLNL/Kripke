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
}
Subdomain::~Subdomain(){
  delete psi;
  delete rhs;
  delete sigt;
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
 * Computes index sets for each octant, and each tile (experimental).
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


/**
 * Grid_Data constructor.
 * Currently, the spatial grid is calculated so that cells are a uniform
 * length = (xmax - xmin) / nx
 * in each spatial direction.
 *
*/
Grid_Data::Grid_Data(Input_Variables *input_vars)
{
  /* Set the processor grid dimensions */
  int R = (input_vars->npx)*(input_vars->npy)*(input_vars->npz);

  /* Check size of PQR_group is the same as MPI_COMM_WORLD */
  int size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if(R != size){
    int myid;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    if(myid == 0){
      printf("ERROR: Incorrect number of MPI tasks. Need %d MPI tasks.", R);
    }
    MPI_Abort(MPI_COMM_WORLD, 1);
  }


  int npx = input_vars->npx;
  int npy = input_vars->npy;
  int npz = input_vars->npz;

  /* Compute the local coordinates in the processor decomposition */
  int myid;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  int isub_ref = myid % npx;
  int jsub_ref = ((myid - isub_ref) / npx) % npy;
  int ksub_ref = (myid - isub_ref - npx*jsub_ref) / (npx * npy);

  /* Compute the processor neighbor array given the above ordering */
  // -1 is used to denote a boundary
  mynbr[0][0] = (isub_ref == 0)     ? -1 : myid-1;
  mynbr[0][1] = (isub_ref == npx-1) ? -1 : myid+1;
  mynbr[1][0] = (jsub_ref == 0)     ? -1 : myid-npx;
  mynbr[1][1] = (jsub_ref == npy-1) ? -1 : myid+npx;
  mynbr[2][0] = (ksub_ref == 0)     ? -1 : myid-npx*npy;
  mynbr[2][1] = (ksub_ref == npz-1) ? -1 : myid+npx*npy;


  // create the kernel object based on nesting
  kernel = createKernel(input_vars->nesting, 3);

  // Create base quadrature set
  InitDirections(this, input_vars->num_dirsets_per_octant * input_vars->num_dirs_per_dirset);

  num_direction_sets = 8*input_vars->num_dirsets_per_octant;
  num_directions_per_set = input_vars->num_dirs_per_dirset;
  num_group_sets = input_vars->num_groupsets;
  num_groups_per_set = input_vars->num_groups_per_groupset;
  num_zone_sets = 1;



  int nx_g = input_vars->nx;
  int ny_g = input_vars->ny;
  int nz_g = input_vars->nz;


  num_legendre = input_vars->legendre_order;
  total_num_moments = (num_legendre+1)*(num_legendre+1);

  int num_subdomains = num_direction_sets*num_group_sets*num_zone_sets;

  Nesting_Order nest = input_vars->nesting;

  // Initialize Subdomains
  subdomains.resize(num_subdomains);
  int group0 = 0;
  for(int gs = 0;gs < num_group_sets;++ gs){
    int dir0 = 0;
    for(int ds = 0;ds < num_direction_sets;++ ds){
      for(int zs = 0;zs < num_zone_sets;++ zs){
        int sdom_id = gs*num_direction_sets*num_zone_sets +
                   ds*num_zone_sets +
                   zs;

        Subdomain &sdom = subdomains[sdom_id];

        // set the set indices
        sdom.idx_group_set = gs;
        sdom.idx_dir_set = ds;
        sdom.idx_zone_set = zs;

        sdom.num_groups = input_vars->num_groups_per_groupset;
        sdom.group0 = group0;

        sdom.num_directions = input_vars->num_dirs_per_dirset;
        sdom.direction0 = dir0;
        sdom.directions = &directions[dir0];

        sdom.deltas[0] = computeGrid(0, npx, nx_g, isub_ref, 0.0, 1.0);
        sdom.deltas[1] = computeGrid(1, npy, ny_g, jsub_ref, 0.0, 1.0);
        sdom.deltas[2] = computeGrid(2, npz, nz_g, ksub_ref, 0.0, 1.0);

        sdom.nzones[0] = sdom.deltas[0].size()-2;
        sdom.nzones[1] = sdom.deltas[1].size()-2;
        sdom.nzones[2] = sdom.deltas[2].size()-2;
        sdom.num_zones = sdom.nzones[0] * sdom.nzones[1] * sdom.nzones[2];

        // allocate storage for the sweep boundary data
        sdom.plane_data[0].resize(sdom.nzones[1] * sdom.nzones[2] * sdom.num_directions * sdom.num_groups);
        sdom.plane_data[1].resize(sdom.nzones[0] * sdom.nzones[2] * sdom.num_directions * sdom.num_groups);
        sdom.plane_data[2].resize(sdom.nzones[0] * sdom.nzones[1] * sdom.num_directions * sdom.num_groups);

        // allocate the storage for solution and source terms
        sdom.psi = new SubTVec(nest, sdom.num_groups, sdom.num_directions, sdom.num_zones);
        sdom.rhs = new SubTVec(nest, sdom.num_groups, sdom.num_directions, sdom.num_zones);
        sdom.sigt = new SubTVec(kernel->nestingSigt(), sdom.num_groups, 1, sdom.num_zones);

        sdom.computeSweepIndexSet();

        // Setup neighbor data
        int dirs[3] = { sdom.directions[0].id, sdom.directions[0].jd, sdom.directions[0].kd};
        for(int dim = 0;dim < 3;++ dim){
          if(dirs[dim] > 0){
            sdom.downwind[dim].mpi_rank = mynbr[dim][1];
            sdom.downwind[dim].subdomain_id = sdom_id;
            sdom.upwind[dim].mpi_rank = mynbr[dim][0];
            sdom.upwind[dim].subdomain_id = sdom_id;
          }
          else{
            sdom.downwind[dim].mpi_rank = mynbr[dim][0];
            sdom.downwind[dim].subdomain_id = sdom_id;
            sdom.upwind[dim].mpi_rank = mynbr[dim][1];
            sdom.upwind[dim].subdomain_id = sdom_id;
          }
        }

      }
      dir0 += input_vars->num_dirs_per_dirset;
    }

    group0 += input_vars->num_groups_per_groupset;
  }

  /* Set ncalls */
  niter = input_vars->niter;

  // setup cross-sections
  sigma_tot.resize(num_group_sets*num_groups_per_set, 0.0);

  // Allocate moments variables
  int total_dirs = num_directions_per_set * num_direction_sets;
  int num_groups = num_groups_per_set * num_group_sets;

  phi = new SubTVec(nest, num_groups, total_num_moments, subdomains[0].num_zones);
  phi_out = new SubTVec(nest, num_groups, total_num_moments, subdomains[0].num_zones);

  ell = new SubTVec(kernel->nestingEll(), total_num_moments, total_dirs, 1);
  ell_plus = new SubTVec(kernel->nestingEllPlus(), total_num_moments, total_dirs, 1);

}

Grid_Data::~Grid_Data(){
  delete kernel;
  delete phi;
  delete phi_out;
  delete ell;
  delete ell_plus;
}

/**
 * Randomizes all variables and matrices for testing suite.
 */
void Grid_Data::randomizeData(void){
  for(int i = 0;i < sigma_tot.size();++i){
    sigma_tot[i] = drand48();
  }

  for(int i = 0;i < directions.size();++i){
    directions[i].xcos = drand48();
    directions[i].ycos = drand48();
    directions[i].zcos = drand48();
  }


  for(int s = 0;s < subdomains.size();++ s){
    subdomains[s].randomizeData();
  }

  phi->randomizeData();
  phi_out->randomizeData();
  ell->randomizeData();
  ell_plus->randomizeData();
}

/**
 * Copies all variables and matrices for testing suite.
 * Correctly copies data from one nesting to another.
 */
void Grid_Data::copy(Grid_Data const &b){
  sigma_tot = b.sigma_tot;
  directions = b.directions;

  subdomains.resize(b.subdomains.size());
  for(int s = 0;s < subdomains.size();++ s){
    subdomains[s].copy(b.subdomains[s]);
  }
  phi->copy(*b.phi);
  phi_out->copy(*b.phi_out);
  ell->copy(*b.ell);
  ell_plus->copy(*b.ell_plus);
}

/**
 * Compares all variables and matrices for testing suite.
 * Correctly compares data from one nesting to another.
 */
bool Grid_Data::compare(Grid_Data const &b, double tol, bool verbose){
  bool is_diff = false;




  for(int i = 0;i < directions.size();++i){
    std::stringstream dirname;
    dirname << "directions[" << i << "]";

    is_diff |= compareScalar(dirname.str()+".xcos",
        directions[i].xcos, b.directions[i].xcos, tol, verbose);

    is_diff |= compareScalar(dirname.str()+".ycos",
        directions[i].ycos, b.directions[i].ycos, tol, verbose);

    is_diff |= compareScalar(dirname.str()+".zcos",
        directions[i].zcos, b.directions[i].zcos, tol, verbose);
  }

  for(int s = 0;s < subdomains.size();++ s){
    is_diff |= subdomains[s].compare(
        b.subdomains[s], tol, verbose);

  }
  is_diff |= compareVector("sigma_tot", sigma_tot, b.sigma_tot, tol, verbose);
  is_diff |= phi->compare("phi", *b.phi, tol, verbose);
  is_diff |= phi_out->compare("phi_out", *b.phi_out, tol, verbose);
  is_diff |= ell->compare("ell", *b.ell, tol, verbose);
  is_diff |= ell_plus->compare("ell_plus", *b.ell_plus, tol, verbose);

  return is_diff;
}


/**
 * Computes the current MPI task's grid given the size of the mesh, and
 * the current tasks index in that dimension (isub_ref).
 */
std::vector<double> Grid_Data::computeGrid(int dim, int npx, int nx_g, int isub_ref, double xmin, double xmax){
 /* Calculate unit roundoff and load into grid_data */
  double eps = 1e-32;
  double thsnd_eps = 1000.e0*(eps);

  // Compute subset of global zone indices
  int nx_l = nx_g / npx;
  int rem = nx_g % npx;
  int ilower, iupper;
  if(rem != 0){
    if(isub_ref < rem){
      nx_l++;
      ilower = isub_ref * nx_l;
    }
    else {
      ilower = rem + isub_ref * nx_l;
    }
  }
  else {
    ilower = isub_ref * nx_l;
  }

  iupper = ilower + nx_l - 1;

  // allocate grid deltas
  std::vector<double> deltas(nx_l+2);

  // Compute the spatial grid
  double dx = (xmax - xmin) / nx_g;
  double coord_lo = xmin + (ilower) * dx;
  double coord_hi = xmin + (iupper+1) * dx;
  for(int i = 0; i < nx_l+2; i++){
    deltas[i] = dx;
  }
  if(std::abs(coord_lo - xmin) <= thsnd_eps*std::abs(xmin)){
    deltas[0] = 0.0;
  }
  if(std::abs(coord_hi - xmax) <= thsnd_eps*std::abs(xmax)){
    deltas[nx_l+1] = 0.0;
  }

  return deltas;
}


