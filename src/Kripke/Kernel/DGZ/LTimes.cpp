#include "../Param.h"

static void dgz_nmd_LTimes_1d(dgz_nmd_Param &p) {
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
      for(int group = 0; group < num_groups; ++group){
        double * __restrict__ psi_d_g = psi_d[group];
        double * __restrict__ phi_n_g = phi_n[group];
        for(int z = 0; z < num_zones; z++){
          double psi_d_g_z = psi_d_g[z];
          phi_n_g[z] += ell_n_m_d * psi_d_g_z;
        }
      }
    }
  }
}

#include <stdio.h>
static void dgz_nmd_LTimes_2d(dgz_nmd_Param &p) {
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
        for(int group = 0; group < num_groups; ++group){
          double *  psi_d_g = psi_d[group];
          double *  phi_n_m_g = phi_n_m[group];
          for(int z = 0; z < num_zones; z++){
            double psi_d_g_z = psi_d_g[z];
            phi_n_m_g[z] += ell_n_m_d * psi_d_g_z;
          }
        }
      }
    }
  }
}


static void dgz_nmd_LTimes_3d(dgz_nmd_Param &p) {
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
        for(int group = 0; group < num_groups; ++group){
          double *  psi_d_g = psi_d[group];
          double *  phi_n_m_g = phi_n_m[group];
          for(int z = 0; z < num_zones; z++){
            double psi_d_g_z = psi_d_g[z];
            phi_n_m_g[z] += ell_n_m_d * psi_d_g_z;
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
void dgz_nmd_Param::LTimes(void) {

  if(geometry_type == 1){
    dgz_nmd_LTimes_1d(*this);
  }
  else if(geometry_type == 2){
    dgz_nmd_LTimes_2d(*this);
  }
  else {
    dgz_nmd_LTimes_3d(*this);
  }
}
