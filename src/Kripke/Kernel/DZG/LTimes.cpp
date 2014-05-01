#include "../Param.h"
#include <algorithm>

static void dzg_nmd_LTimes_1d(dzg_nmd_Param &p) {
  // Grab parameters
  double ***psi = p.psi.data;
  double ***phi = p.phi.data;
  double ***ell = p.ell.data;

  int num_groups = p.num_groups;
  int num_local_directions = p.num_directions;
  int num_zones = p.num_zones;

  int num_moments = p.num_moments;

  /* 1D Spherical or Slab Geometry */
  p.phi.clear(0.0);
  for(int n = 0; n < num_moments; n++){
    double **phi_n = phi[n];
    for(int d = 0; d < num_local_directions; d++){
      double **psi_d = psi[d];
      double ell_n_m_d = ell[n][0][d];

#ifdef _OPENMP
#pragma omp parallel for
#endif
      for(int z = 0; z < num_zones; z++){
        double * __restrict__ psi_d_z = psi_d[z];
        double * __restrict__ phi_n_z = phi_n[z];
        for(int group = 0; group < num_groups; ++group){
          double psi_d_g_z = psi_d_z[group];
          phi_n_z[group] += ell_n_m_d * psi_d_g_z;
        }
      }
    }
  }
}

#include <stdio.h>
static void dzg_nmd_LTimes_2d(dzg_nmd_Param &p) {
  // Grab parameters
  double ***psi = p.psi.data;
  double ***phi = p.phi.data;
  double ***ell = p.ell.data;

  int num_groups = p.num_groups;
  int num_local_directions = p.num_directions;
  int num_zones = p.num_zones;

  int num_moments = p.num_moments;

  /* 2D rho-z Cylindrical Geometry */
  p.phi.clear(0.0);
  for(int n = 0; n < num_moments; n++){
    double ***phi_n = phi + (n*(n+1)/2);
    double **ell_n = ell[n];
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(int m = 0; m <= n; m++){
      double **phi_n_m = phi_n[m];
      double *ell_n_m = ell_n[m];
      for(int d = 0; d < num_local_directions; d++){
        double **psi_d = psi[d];
        double ell_n_m_d = ell_n_m[d];
        for(int z = 0; z < num_zones; z++){
          double * __restrict__ psi_d_z = psi_d[z];
          double * __restrict__ phi_n_z = phi_n_m[z];
          for(int group = 0; group < num_groups; ++group){
            double psi_d_g_z = psi_d_z[group];
            phi_n_z[group] += ell_n_m_d * psi_d_g_z;
          }
        }
      }
    }
  }
}


static void dzg_nmd_LTimes_3d(dzg_nmd_Param &p) {
  // Grab parameters
  double ***psi = p.psi.data;
  double ***phi = p.phi.data;
  double ***ell = p.ell.data;

  int num_groups = p.num_groups;
  int num_local_directions = p.num_directions;
  int num_zones = p.num_zones;

  int num_moments = p.num_moments;

  /* 3D Cartesian Geometry */
  p.phi.clear(0.0);
  for(int n = 0; n < num_moments; n++){
    double ***phi_n = phi + n*n;
    double **ell_n = ell[n];
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(int m = -n; m <= n; m++){
      double **phi_n_m = phi_n[m+n];
      double *ell_n_m = ell_n[m+n];
      for(int d = 0; d < num_local_directions; d++){
        double **psi_d = psi[d];
        double ell_n_m_d = ell_n_m[d];
        for(int z = 0; z < num_zones; z++){
          double * __restrict__ psi_d_z = psi_d[z];
          double * __restrict__ phi_n_z = phi_n_m[z];
          for(int group = 0; group < num_groups; ++group){
            double psi_d_g_z = psi_d_z[group];
            phi_n_z[group] += ell_n_m_d * psi_d_g_z;
          }
        }
      }
    }
  }
}


/**
 * LTimes routine with data nesting of:
 *    Psi[group][direction][zone]
 *    Phi[moment][zone]
 */
void dzg_nmd_Param::LTimes(void) {

  if(geometry_type == 1){
    dzg_nmd_LTimes_1d(*this);
  }
  else if(geometry_type == 2){
    dzg_nmd_LTimes_2d(*this);
  }
  else {
    dzg_nmd_LTimes_3d(*this);
  }
}
