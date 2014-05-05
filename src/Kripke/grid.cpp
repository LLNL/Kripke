/*--------------------------------------------------------------------------
 * Utility functions for the Grid_Data structure.
 *--------------------------------------------------------------------------*/

#include <Kripke/grid.h>

#include <Kripke/comm.h>
#include <Kripke/input_variables.h>

#include <cmath>

Group_Dir_Set::Group_Dir_Set() :
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
Group_Dir_Set::~Group_Dir_Set(){
  delete psi;
  delete rhs;
  delete sigt;
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
  
  computeGrid(0, npx, nx_g, isub_ref, input_vars->xmin, input_vars->xmax);
  computeGrid(1, npy, ny_g, jsub_ref, input_vars->ymin, input_vars->ymax);
  computeGrid(2, npz, nz_g, ksub_ref, input_vars->zmin, input_vars->zmax);
  num_zones = nzones[0]*nzones[1]*nzones[2];
  tmp_sigma_tot.resize(num_zones, 0.0);
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
