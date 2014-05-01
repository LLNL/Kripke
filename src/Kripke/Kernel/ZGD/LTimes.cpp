#include "../Param.h"

static void zgd_nmd_LTimes_1d(zgd_nmd_Param &p) {
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
  for(int z = 0; z < num_zones; z++){
    double **psi_z = psi[z];
    double **phi_z = phi[z];

    for(int group = 0; group < num_groups; ++group){
      double * __restrict__ psi_z_g = psi_z[group];
      double * __restrict__ phi_z_g = phi_z[group];

      for(int n = 0; n < num_moments; n++){
        double * __restrict__ ell_n_m = ell[n][0]; // m == 0

        double phi_z_g_n = 0.;
        for(int d = 0; d < num_local_directions; d++){
          double ell_n_m_d = ell_n_m[d];
          double psi_z_g_d = psi_z_g[d];
          phi_z_g_n += ell_n_m_d * psi_z_g_d;
        }

        phi_z_g[n] = phi_z_g_n;
      }
    }
  }

}


static void zgd_nmd_LTimes_2d(zgd_nmd_Param &p) {
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
  for(int z = 0; z < num_zones; z++){
    double **psi_z = psi[z];
    double **phi_z = phi[z];

    for(int group = 0; group < num_groups; ++group){
      double * __restrict__ psi_z_g = psi_z[group];
      double * __restrict__ phi_z_g = phi_z[group];

      for(int n = 0; n < num_moments; n++){
        double * __restrict__ phi_z_g_n = phi_z_g +(n*(n+1)/2);
        double **ell_n = ell[n];

        for(int m = 0; m <= n; m++){
          double * __restrict__ ell_n_m = ell[n][m];

          double phi_z_g_n_m = 0.;
          for(int d = 0; d < num_local_directions; d++){
            double ell_n_m_d = ell_n_m[d];
            double psi_z_g_d = psi_z_g[d];
            phi_z_g_n_m += ell_n_m_d * psi_z_g_d;
          }

          phi_z_g_n[m] = phi_z_g_n_m;
        }
      }
    }
  }
}


static void zgd_nmd_LTimes_3d(zgd_nmd_Param &p) {
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
  for(int z = 0; z < num_zones; z++){
    double **psi_z = psi[z];
    double **phi_z = phi[z];

    for(int group = 0; group < num_groups; ++group){
      double * __restrict__ psi_z_g = psi_z[group];
      double * __restrict__ phi_z_g = phi_z[group];

      for(int n = 0; n < num_moments; n++){
        double * __restrict__ phi_z_g_n = phi_z_g + n*n + n;
        double **ell_n = ell[n];

        for(int m = -n; m <= n; m++){
          double * __restrict__ ell_n_m = ell[n][m+n];

          double phi_z_g_n_m = 0.;
          for(int d = 0; d < num_local_directions; d++){
            double ell_n_m_d = ell_n_m[d];
            double psi_z_g_d = psi_z_g[d];
            phi_z_g_n_m += ell_n_m_d * psi_z_g_d;
          }

          phi_z_g_n[m] = phi_z_g_n_m;
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
void zgd_nmd_Param::LTimes(void) {

  if(geometry_type == 1){
    zgd_nmd_LTimes_1d(*this);
  }
  else if(geometry_type == 2){
    zgd_nmd_LTimes_2d(*this);
  }
  else {
    zgd_nmd_LTimes_3d(*this);
  }
}
