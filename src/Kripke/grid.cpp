/*--------------------------------------------------------------------------
 * Utility functions for the Grid_Data structure.
 *--------------------------------------------------------------------------*/

#include "transport_headers.h"
#include <cmath>
#include <float.h>

/*--------------------------------------------------------------------------
 * GenGrid : Creates a new Grid_Data structure and allocates
 *                 memory for its data.
 *
 * Currently, the spatial grid is calculated so that cells are a uniform
 * length = (xmax - xmin) / nx
 * in each spatial direction.
 *
 *--------------------------------------------------------------------------*/

Grid_Data::Grid_Data(int npx, int npy, int npz, int num_directions_per_octant,
                   double xmin, double xmax, int nx_g,
                   double ymin, double ymax, int ny_g,
                   double zmin, double zmax, int nz_g,
                   MPI_Comm comm)
/*--------------------------------------------------------------------------
 * npx            : The number of processors in the x direction
 * npy            : The number of processors in the y direction
 * npz            : The number of processors in the z direction
 * num_directions_per_octant : The number of directions in each octant
 * xmin           : The miniumum spatial value for the x direction
 * xmax           : The maximum spatial value for the x direction
 * nx_g           : The number of spatial zones to use in the x direction
 * ymin           : The miniumum spatial value for the y direction
 * ymax           : The maximum spatial value for the y direction
 * ny_g           : The number of spatial zones to use in the y direction
 * zmin           : The miniumum spatial value for the z direction
 * zmax           : The maximum spatial value for the z direction
 * nz_g           : The number of spatial zones to use in the z direction
 * comm           : Spatial communicator
 *--------------------------------------------------------------------------*/
{
  std::vector<double> deltas_x, deltas_y, deltas_z;
  double dx, dy, dz;

  int myid;
  int isub_ref, jsub_ref, ksub_ref;
  int nx_l, ny_l, nz_l;
  int rem;
  int i;
  int lmin, lmax;

  /* Calculate unit roundoff and load into grid_data */
  eps = DBL_EPSILON;
  double thsnd_eps = 1000.e0*(eps);

  /* Compute the local coordinates in the processor decomposition */
  MPI_Comm_rank(comm, &myid);

  nprocs[0] = npx;
  nprocs[1] = npy;
  nprocs[2] = npz;

  isub_ref = myid % npx;
  jsub_ref = ((myid - isub_ref) / npx) % npy;
  ksub_ref = (myid - isub_ref - npx*jsub_ref) / (npx * npy);

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

  /* Compute the x-direction local extents */
  nx_l = nx_g / npx;
  rem = nx_g % npx;
  if(rem != 0){
    if(isub_ref < rem){
      nx_l++;
      ilower[0] = isub_ref * nx_l;
    }
    else {
      ilower[0] = rem + isub_ref * nx_l;
    }
  }
  else {
    ilower[0] = isub_ref * nx_l;
  }

  iupper[0] = ilower[0] + nx_l - 1;

  /* Compute the y-direction local extents */
  ny_l = ny_g / npy;
  rem = ny_g % npy;
  if(rem != 0){
    if(jsub_ref < rem){
      ny_l++;
      ilower[1] = jsub_ref * ny_l;
    }
    else {
      ilower[1] = rem + jsub_ref * ny_l;
    }
  }
  else {
    ilower[1] = jsub_ref * ny_l;
  }

  iupper[1] = ilower[1] + ny_l - 1;

  /* Compute the z-direction local extents */
  nz_l = nz_g / npz;
  rem = nz_g % npz;
  if(rem != 0){
    if(ksub_ref < rem){
      nz_l++;
      ilower[2] = ksub_ref * nz_l;
    }
    else {
      ilower[2] = rem + ksub_ref * nz_l;
    }
  }
  else {
    ilower[2] = ksub_ref * nz_l;
  }

  iupper[2] = ilower[2] + nz_l - 1;

  /* Allocate the spatial grid */
  deltas_x.resize(nx_l+2);
  deltas_y.resize(ny_l+2);
  deltas_z.resize(nz_l+2);

  /* Compute the x-direction spatial grid */
  double coord_lo[3] = {0, 0, 0};
  double coord_hi[3] = {0, 0, 0};
  dx = (xmax - xmin) / nx_g;
  coord_lo[0] = xmin + (ilower[0]) * dx;
  coord_hi[0] = xmin + (iupper[0]+1) * dx;
  for(i = 0; i < nx_l+2; i++){
    deltas_x[i] = dx;
  }
  if(std::abs(coord_lo[0] - xmin) <= thsnd_eps*std::abs(xmin)){
    deltas_x[0] = 0.0;
  }
  if(std::abs(coord_hi[0] - xmax) <= thsnd_eps*std::abs(xmax)){
    deltas_x[nx_l+1] = 0.0;
  }

  /* Compute the y-direction spatial grid */
  dy = (ymax - ymin) / ny_g;
  coord_lo[1] = ymin + (ilower[1]) * dy;
  coord_hi[1] = ymin + (iupper[1]+1) * dy;
  for(i = 0; i < ny_l+2; i++){
    deltas_y[i] = dy;
  }
  if(std::abs(coord_lo[1] - ymin) <= thsnd_eps*std::abs(ymin)){
    deltas_y[0] = 0.0;
  }
  if(std::abs(coord_hi[1] - ymax) <= thsnd_eps*std::abs(ymax)){
    deltas_y[ny_l+1] = 0.0;
  }

  /* Compute the z-direction spatial grid */
  dz = (zmax - zmin) / nz_g;
  coord_lo[2] = zmin + (ilower[2]) * dz;
  coord_hi[2] = zmin + (iupper[2]+1) * dz;
  for(i = 0; i < nz_l+2; i++){
    deltas_z[i] = dz;
  }
  if(std::abs(coord_lo[2] - zmin) <= thsnd_eps*std::abs(zmin)){
    deltas_z[0] = 0.0;
  }
  if(std::abs(coord_hi[2] - zmax) <= thsnd_eps*std::abs(zmax)){
    deltas_z[nz_l+1] = 0.0;
  }

  /* Compute the volumes */
  int num_zones = nx_l*ny_l*nz_l;
  volume.resize(num_zones);
  for(int k = 0; k < nz_l; k++){
    for(int j = 0; j < ny_l; j++){
      for(int i = 0; i < nx_l; i++){
        int ijk = i + nx_l*j + nx_l*ny_l*k;
        volume[ijk] = deltas_x[i+1]*deltas_y[j+1]*deltas_z[k+1];
      }
    }
  }

  nzones[0] = nx_l;
  nzones[1] = ny_l;
  nzones[2] = nz_l;

  tmp_sigma_tot.resize(num_zones, 0.0);

  InitDirections(this, num_directions_per_octant);

  global_num_zones      = (double)nx_g * (double)ny_g * (double)nz_g;
  deltas[0]             = deltas_x;
  deltas[1]             = deltas_y;
  deltas[2]             = deltas_z;

  int num_directions = directions.size();
  int i_plane_zones = ny_l * nz_l;
  int j_plane_zones = nx_l * nz_l;
  int k_plane_zones = nx_l * ny_l;
  psi_i_plane.resize(num_directions);
  psi_j_plane.resize(num_directions);
  psi_k_plane.resize(num_directions);
  for(int d=0; d<num_directions; d++){
    psi_i_plane[d].resize(i_plane_zones);
    psi_j_plane[d].resize(j_plane_zones);
    psi_k_plane[d].resize(k_plane_zones);
  }

}

