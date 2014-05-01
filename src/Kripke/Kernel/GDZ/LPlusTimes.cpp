#include "../Param.h"

void gdz_nmd_LPlusTimes_1d(gdz_nmd_Param &p) {
  // Grab parameters
  double ***psi = p.psi.data;
  double ***phi = p.phi.data;
  double ***ell_plus = p.ell_plus.data;

  int num_groups = p.num_groups;
  int num_local_directions = p.num_directions;
  int num_zones = p.num_zones;

  int num_moments = p.num_moments;

  /* 1D Spherical or Slab Geometry */
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(int group = 0; group < num_groups; ++group){
    double **phi_g = phi[group];
    double **psi_g = psi[group];

    for(int d=0; d < num_local_directions; d++){
      double **ell_plus_d = ell_plus[d];
      double * __restrict__ psi_g_d = psi_g[d];
      for(int i=0; i < num_zones; i++){
        psi_g_d[i] = 0.0;
      }

      for(int n=0; n< num_moments; n++){
        double * __restrict__ phi_g_n = phi_g[n];
        double *ell_plus_d_n = ell_plus_d[n];
        double ell_plus_d_n_m = ell_plus_d_n[0]; // m = 0;
        for(int i=0; i < num_zones; i++){
          psi_g_d[i] += ell_plus_d_n_m * phi_g_n[i];
        }
      }
    }
  }
}


void gdz_nmd_LPlusTimes_2d(gdz_nmd_Param &p) {
  // Grab parameters
  double ***psi = p.psi.data;
  double ***phi = p.phi.data;
  double ***ell_plus = p.ell_plus.data;

  int num_groups = p.num_groups;
  int num_local_directions = p.num_directions;
  int num_zones = p.num_zones;

  int num_moments = p.num_moments;

  /* 2D rho-z Cylindrical Geometry */
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(int group = 0; group < num_groups; ++group){
    double **phi_g = phi[group];
    double **psi_g = psi[group];

    for(int d=0; d < num_local_directions; d++){
      double **ell_plus_d = ell_plus[d];
      double * __restrict__ psi_g_d = psi_g[d];
      for(int i=0; i < num_zones; i++){
        psi_g_d[i] = 0.0;
      }

      for(int n=0; n< num_moments; n++){
        double **phi_g_n = phi_g + n*(n+1)/2;
        double *ell_plus_d_n = ell_plus_d[n];

        for(int m = 0; m <= n; m++){
          double ell_plus_d_n_m = ell_plus_d_n[m];
          double * __restrict__ phi_g_nm = phi_g_n[m];

          for(int i=0; i < num_zones; i++){
            psi_g_d[i] += ell_plus_d_n_m * phi_g_nm[i];
          }
        }
      }
    }
  }
}


void gdz_nmd_LPlusTimes_3d(gdz_nmd_Param &p) {
  // Grab parameters
  double ***psi = p.psi.data;
  double ***phi = p.phi.data;
  double ***ell_plus = p.ell_plus.data;

  int num_groups = p.num_groups;
  int num_local_directions = p.num_directions;
  int num_zones = p.num_zones;

  int num_moments = p.num_moments;

  /* 3D Cartesian Geometry */
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(int group = 0; group < num_groups; ++group){
    double **phi_g = phi[group];
    double **psi_g = psi[group];

    for(int d=0; d < num_local_directions; d++){
      double **ell_plus_d = ell_plus[d];
      double * __restrict__ psi_g_d = psi_g[d];
      for(int i=0; i < num_zones; i++){
        psi_g_d[i] = 0.0;
      }

      for(int n=0; n< num_moments; n++){
        double **phi_g_n = phi_g + n*n;
        double *ell_plus_d_n = ell_plus_d[n];

        for(int m = -n; m <= n; m++){
          double ell_plus_d_n_m = ell_plus_d_n[m+n];
          double * __restrict__ phi_g_nm = phi_g_n[m+n];

          for(int i=0; i < num_zones; i++){
            psi_g_d[i] += ell_plus_d_n_m * phi_g_nm[i];
          }
        }
      }
    }
  }
}


/**
 * LPlusTimes routine with data nesting of:
 *    Psi[group][direction][zone]
 *    Phi[moment][zone]
 */
void gdz_nmd_Param::LPlusTimes(void) {

  if(geometry_type == 1){
    gdz_nmd_LPlusTimes_1d(*this);
  }
  else if(geometry_type == 2){
    gdz_nmd_LPlusTimes_2d(*this);
  }
  else {
    gdz_nmd_LPlusTimes_3d(*this);
  }
}
