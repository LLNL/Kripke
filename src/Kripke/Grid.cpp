/*--------------------------------------------------------------------------
 * Utility functions for the Grid_Data structure.
 *--------------------------------------------------------------------------*/

#include <Kripke/Grid.h>
#include <Kripke/SubTVec.h>
#include <Kripke/LMat.h>
#include <Kripke/Comm.h>
#include <Kripke/Input_Variables.h>

#include <cmath>

Group_Dir_Set::Group_Dir_Set() :
  num_groups(0),
  num_directions(0),
  group0(0),
  direction0(0),
  directions(NULL),
  psi(NULL),
  rhs(NULL),
  sigt(NULL),
  psi_lf(NULL),
  psi_fr(NULL),
  psi_bo(NULL)
{
}
Group_Dir_Set::~Group_Dir_Set(){
  delete psi;
  delete rhs;
  delete sigt;

  delete psi_lf;
  delete psi_fr;
  delete psi_bo;
}


void Group_Dir_Set::allocate(Grid_Data *grid_data, Nesting_Order nest){
  delete psi;
  psi = new SubTVec(nest,
      num_groups, num_directions, grid_data->num_zones);

  delete rhs;
  rhs = new SubTVec(nest,
      num_groups, num_directions, grid_data->num_zones);

  // allocate sigt  1xGxZ if groups come before zones
  delete sigt;
  if(nest == NEST_GDZ || nest ==  NEST_DGZ || nest == NEST_GZD){
    sigt = new SubTVec(NEST_DGZ,
      num_groups, 1, grid_data->num_zones);
  }
  // otherwise, 1xZxG
  else{
    sigt = new SubTVec(NEST_DZG,
      num_groups, 1, grid_data->num_zones);
  }

  // Allocate sweep boundary data
  int local_imax = grid_data->nzones[0];
  int local_jmax = grid_data->nzones[1];
  int local_kmax = grid_data->nzones[2];
  psi_lf = new SubTVec(nest, num_groups, num_directions,
                    (local_imax+1)*local_jmax*local_kmax);
  psi_fr = new SubTVec(nest, num_groups, num_directions,
                    local_imax*(local_jmax+1)*local_kmax);
  psi_bo = new SubTVec(nest, num_groups, num_directions,
                    local_imax*local_jmax*(local_kmax+1));
}


/*--------------------------------------------------------------------------
 * GenGrid : Creates a new Grid_Data structure and allocates
 *                 memory for its data.
 *
 * Currently, the spatial grid is calculated so that cells are a uniform
 * length = (xmax - xmin) / nx
 * in each spatial direction.
 *
 *--------------------------------------------------------------------------*/

Grid_Data::Grid_Data(Input_Variables *input_vars, Directions *directions)
{
  int npx = input_vars->npx;
  int npy = input_vars->npy;
  int npz = input_vars->npz;
  int nx_g = input_vars->nx;
  int ny_g = input_vars->ny;
  int nz_g = input_vars->nz;

  /* Compute the local coordinates in the processor decomposition */
  int myid;
  MPI_Comm_rank(GetRGroup(), &myid);

  int isub_ref = myid % npx;
  int jsub_ref = ((myid - isub_ref) / npx) % npy;
  int ksub_ref = (myid - isub_ref - npx*jsub_ref) / (npx * npy);

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
  
  computeGrid(0, npx, nx_g, isub_ref, 0.0, 1.0);
  computeGrid(1, npy, ny_g, jsub_ref, 0.0, 1.0);
  computeGrid(2, npz, nz_g, ksub_ref, 0.0, 1.0);
  num_zones = nzones[0]*nzones[1]*nzones[2];

  num_moments = 4;

  sig_s.resize(num_zones, 0.0);
}

Grid_Data::~Grid_Data(){
  delete phi;
  delete phi_out;
  delete ell;
  delete ell_plus;
}


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
