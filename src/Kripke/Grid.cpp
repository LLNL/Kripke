#include <Kripke/Grid.h>
#include <Kripke/SubTVec.h>
#include <Kripke/Comm.h>
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
  rhs(NULL)
{
}
Subdomain::~Subdomain(){
  delete psi;
  delete rhs;
}


/**
 * Randomizes data for a set.
 */
void Subdomain::randomizeData(void){
  psi->randomizeData();
  rhs->randomizeData();
}

/**
 * Copies two sets, allowing for different nestings.
 */
void Subdomain::copy(Subdomain const &b){
  psi->copy(*b.psi);
  rhs->copy(*b.rhs);
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

  return is_diff;
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
  int R = (input_vars->npx)*(input_vars->npy)*(input_vars->npz);;
  /* Check size of PQR_group is the same as MPI_COMM_WORLD */
  int size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if(R != size){
    int myid;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    if(myid == 0){
      printf("ERROR: Incorrect number of MPI tasks. Need %d MPI tasks.", R);
    }
    error_exit(1);
  }

  /* Compute the local coordinates in the processor decomposition */
  int npx = input_vars->npx;
  int npy = input_vars->npy;
  int npz = input_vars->npz;

  int myid;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  int isub_ref = myid % npx;
  int jsub_ref = ((myid - isub_ref) / npx) % npy;
  int ksub_ref = (myid - isub_ref - npx*jsub_ref) / (npx * npy);


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

  computeGrid(0, npx, nx_g, isub_ref, 0.0, 1.0);
  computeGrid(1, npy, ny_g, jsub_ref, 0.0, 1.0);
  computeGrid(2, npz, nz_g, ksub_ref, 0.0, 1.0);

  nzones[0] = deltas[0].size()-2;
  nzones[1] = deltas[1].size()-2;
  nzones[2] = deltas[2].size()-2;
  num_zones= nzones[0] * nzones[1] * nzones[2];

  num_legendre = input_vars->legendre_order;
  total_num_moments = (num_legendre+1)*(num_legendre+1);

  int num_subdomains = num_direction_sets*num_group_sets*num_zone_sets;

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

        sdom.num_zones = num_zones;

        // allocate the storage for solution and source terms
        Nesting_Order nest = kernel->nestingPsi();
        sdom.psi = new SubTVec(nest, sdom.num_groups, sdom.num_directions, sdom.num_zones);
        sdom.rhs = new SubTVec(nest, sdom.num_groups, sdom.num_directions, sdom.num_zones);

      }
      dir0 += input_vars->num_dirs_per_dirset;
    }

    group0 += input_vars->num_groups_per_groupset;
  }

  /* Set ncalls */
  niter = input_vars->niter;

  // setup cross-sections
  sigma_tot.resize(num_group_sets*num_groups_per_set, 0.0);

  kernel->allocateStorage(this);

  /* Create buffer info for sweeping if using Diamond-Difference */
  CreateBufferInfo(this);





  /* Compute the processor neighbor array assuming a lexigraphic ordering */
  if(isub_ref == 0){
    mynbr[0][0] = -1;
  }
  else {
    mynbr[0][0] = myid - 1;
  }

  if(isub_ref == npx-1){
    mynbr[0][1] = -1;
  }
  else {
    mynbr[0][1] = myid + 1;
  }

  if(jsub_ref == 0){
    mynbr[1][0] = -1;
  }
  else {
    mynbr[1][0] = myid - npx;
  }

  if(jsub_ref == npy-1){
    mynbr[1][1] = -1;
  }
  else {
    mynbr[1][1] = myid + npx;
  }

  if(ksub_ref == 0){
    mynbr[2][0] = -1;
  }
  else {
    mynbr[2][0] = myid - npx * npy;
  }

  if(ksub_ref == npz-1){
    mynbr[2][1] = -1;
  }
  else {
    mynbr[2][1] = myid + npx * npy;
  }





  computeSweepIndexSets(input_vars->block_size);
}

Grid_Data::~Grid_Data(){
  delete comm;
  delete kernel;
  delete sigt;
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

  for(int d = 0;d < 3;++ d){
    for(int i = 0;i < deltas[d].size();++ i){
      deltas[d][i] = drand48();
    }
  }

  for(int s = 0;s < subdomains.size();++ s){
    subdomains[s].randomizeData();
  }

  sigt->randomizeData();
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

  for(int d = 0;d < 3;++ d){
    deltas[d] = b.deltas[d];
  }

  subdomains.resize(b.subdomains.size());
  for(int s = 0;s < subdomains.size();++ s){
    subdomains[s].copy(b.subdomains[s]);
  }
  sigt->copy(*b.sigt);
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


  is_diff |= compareVector("sigma_tot", sigma_tot, b.sigma_tot, tol, verbose);

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

  is_diff |= compareVector("deltas[0]", deltas[0], b.deltas[0], tol, verbose);
  is_diff |= compareVector("deltas[1]", deltas[1], b.deltas[1], tol, verbose);
  is_diff |= compareVector("deltas[2]", deltas[2], b.deltas[2], tol, verbose);

  for(int s = 0;s < subdomains.size();++ s){
    is_diff |= subdomains[s].compare(
        b.subdomains[s], tol, verbose);

  }

  is_diff |= sigt->compare("sigt", *b.sigt, tol, verbose);
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
void Grid_Data::computeGrid(int dim, int npx, int nx_g, int isub_ref, double xmin, double xmax){
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
  deltas[dim].resize(nx_l+2);

  // Compute the spatial grid
  double dx = (xmax - xmin) / nx_g;
  double coord_lo = xmin + (ilower) * dx;
  double coord_hi = xmin + (iupper+1) * dx;
  for(int i = 0; i < nx_l+2; i++){
    deltas[dim][i] = dx;
  }
  if(std::abs(coord_lo - xmin) <= thsnd_eps*std::abs(xmin)){
    deltas[dim][0] = 0.0;
  }
  if(std::abs(coord_hi - xmax) <= thsnd_eps*std::abs(xmax)){
    deltas[dim][nx_l+1] = 0.0;
  }

  nzones[dim] = nx_l;
}

/**
 * Computes index sets for each octant, and each tile (experimental).
 * Determines logical indices, and increments for i,j,k based on grid
 * information and quadrature set sweeping direction.
 */
void Grid_Data::computeSweepIndexSets(int block_size){
  octant_extent.resize(8);
  for(int octant = 0;octant < 8;++ octant){

    int id, jd, kd;
    switch(octant){
      case 0: id = 1; jd = 1; kd = 1; break;
      case 1: id = -1; jd = 1; kd = 1; break;
      case 2: id = -1; jd = -1; kd = 1; break;
      case 3: id = 1; jd = -1; kd = 1; break;
      case 4: id = 1; jd = 1; kd = -1; break;
      case 5: id = -1; jd = 1; kd = -1; break;
      case 6: id = -1; jd = -1; kd = -1; break;
      case 7: id = 1; jd = -1; kd = -1; break;
    }

    int istartz, istopz, in, il, ir;

    if(id > 0){
      istartz = 0; istopz = nzones[0]-1; in = 1; il = 0; ir = 1;
    }
    else {
      istartz = nzones[0]-1; istopz = 0; in = -1; il = 1; ir = 0;
    }

    int jstartz, jstopz, jn, jf, jb;
    if(jd > 0){
      jstartz = 0; jstopz = nzones[1]-1; jn = 1; jf = 0; jb = 1;
    }
    else {
      jstartz = nzones[1]-1; jstopz = 0; jn = -1; jf = 1; jb = 0;
    }

    int kstartz, kstopz, kn, kb, kt;
    if(kd > 0){
      kstartz = 0; kstopz = nzones[2]-1; kn =  1; kb = 0; kt = 1;
    }
    else {
      kstartz = nzones[2]-1; kstopz = 0; kn = -1; kb = 1; kt = 0;
    }

    // Define extent block and pattern
    Grid_Sweep_Block &extent = octant_extent[octant];
    extent.start_i = istartz;
    extent.start_j = jstartz;
    extent.start_k = kstartz;
    extent.end_i = istopz + in;
    extent.end_j = jstopz + jn;
    extent.end_k = kstopz + kn;
    extent.inc_i = in;
    extent.inc_j = jn;
    extent.inc_k = kn;
  }
}
