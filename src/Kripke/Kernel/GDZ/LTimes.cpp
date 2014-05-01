#include "../Param.h"

void gdz_nmd_LTimes_1d(gdz_nmd_Param &p) {
  // Grab parameters
  double ***psi = p.psi.data;
  double ***phi = p.phi.data;
  double ***ell = p.ell.data;

  int num_groups = p.num_groups;
  int num_local_directions = p.num_directions;
  int num_zones = p.num_zones;

  int num_moments = p.num_moments;

  /* 1D Spherical or Slab Geometry */
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(int group = 0; group < num_groups; ++group){
    double **psi_zonal = psi[group];
    double **phi_g = phi[group];
    for(int n = 0; n < num_moments; n++){
      double *phi_n = phi_g[n];
      double **ell_n = ell[n];
      for(int i = 0; i < num_zones; i++){
        phi_n[i] = 0.0;
      }
      int m = 0;
      double *ell_n_m = ell_n[m];
      for(int d = 0; d < num_local_directions; d++){
        double ell_n_m_d = ell_n_m[d];
        double *psi_zonal_d = psi_zonal[d];
        for(int i = 0; i < num_zones; i++){
          phi_n[i] += ell_n_m_d * psi_zonal_d[i];
        }
      }
    }
  }
}


void gdz_nmd_LTimes_2d(gdz_nmd_Param &p) {
  // Grab parameters
  double ***psi = p.psi.data;
  double ***phi = p.phi.data;
  double ***ell = p.ell.data;

  int num_groups = p.num_groups;
  int num_local_directions = p.num_directions;
  int num_zones = p.num_zones;

  int num_moments = p.num_moments;

  /* 2D rho-z Cylindrical Geometry */
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(int group = 0; group < num_groups; ++group){
    double **psi_zonal = psi[group];
    double **phi_g = phi[group];
    for(int n = 0; n < num_moments; n++){
      double **ell_n = ell[n];
      for(int m = 0; m <= n; m++){
        double *phi_n_m = phi_g[n*(n+1)/2 + m];
        for(int i = 0; i < num_zones; i++){
          phi_n_m[i] = 0.0;
        }
        for(int d = 0; d < num_local_directions; d++){
          double *psi_zonal_d = psi_zonal[d];
          double ell_n_m_d = ell_n[m][d];
          for(int i = 0; i < num_zones; i++){
            double val = psi_zonal_d[i];
            phi_n_m[i] += ell_n_m_d * val;
          }
        }
      }
    }
  }
}


void gdz_nmd_LTimes_3d(gdz_nmd_Param &p) {
  // Grab parameters
  double ***psi = p.psi.data;
  double ***phi = p.phi.data;
  double ***ell = p.ell.data;

  int num_groups = p.num_groups;
  int num_local_directions = p.num_directions;
  int num_zones = p.num_zones;

  int num_moments = p.num_moments;

  /* 3D Cartesian Geometry */
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(int group = 0; group < num_groups; ++group){
    double **psi_zonal = psi[group];
    double **phi_g = phi[group];
    for(int n = 0; n < num_moments; n++){
      double **ell_n = ell[n];
      double **phi_g_n = phi_g + n*n + n;
      for(int m = -n; m <= n; m++){
        double *__restrict__ phi_g_nm = phi_g_n[m];
        double * __restrict__ ell_n_m = ell_n[m+n];
        for(int i = 0; i < num_zones; i++){
          phi_g_nm[i] = 0.0;
        }
        for(int d = 0; d < num_local_directions; d++){
          double ell_n_m_d = ell_n_m[d];
          double * __restrict__ psi_g_d = psi_zonal[d];
          for(int i = 0; i < num_zones; i++){
            phi_g_nm[i] += ell_n_m_d * psi_g_d[i];
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
void gdz_nmd_Param::LTimes(void) {

  if(geometry_type == 1){
    gdz_nmd_LTimes_1d(*this);
  }
  else if(geometry_type == 2){
    gdz_nmd_LTimes_2d(*this);
  }
  else {
    gdz_nmd_LTimes_3d(*this);
  }
}
