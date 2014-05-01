#include "../Param.h"

static void zgd_nmd_LPlusTimes_1d(zgd_nmd_Param &p) {
  // Grab parameters
  double ***psi = p.psi.data;
  double ***phi = p.phi.data;
  double ***ell_plus = p.ell_plus.data;

  int num_groups = p.num_groups;
  int num_local_directions = p.num_directions;
  int num_zones = p.num_zones;

  int num_moments = p.num_moments;

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(int z = 0; z < num_zones; z++){
    double **psi_z = psi[z];
    double **phi_z = phi[z];

    for(int group = 0; group < num_groups; ++group){
      double * __restrict__ psi_z_g = psi_z[group];
      double * __restrict__ phi_z_g = phi_z[group];

      for(int d = 0; d < num_local_directions; d++){
        double **ell_plus_d = ell_plus[d];
        double psi_z_g_d = 0.0;

        for(int n = 0; n < num_moments; n++){
          double *ell_plus_d_n = ell_plus_d[n];
          double ell_plus_d_n_m = ell_plus_d_n[0]; // m = 0;
          double phi_z_g_n = phi_z_g[n];

          psi_z_g_d += ell_plus_d_n_m * phi_z_g_n;
        }

        psi_z_g[d] = psi_z_g_d;
      }
    }
  }

}


static void zgd_nmd_LPlusTimes_2d(zgd_nmd_Param &p) {
  // Grab parameters
  double ***psi = p.psi.data;
  double ***phi = p.phi.data;
  double ***ell_plus = p.ell_plus.data;

  int num_groups = p.num_groups;
  int num_local_directions = p.num_directions;
  int num_zones = p.num_zones;

  int num_moments = p.num_moments;

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(int z = 0; z < num_zones; z++){
    double **psi_z = psi[z];
    double **phi_z = phi[z];

    for(int group = 0; group < num_groups; ++group){
      double * __restrict__ psi_z_g = psi_z[group];
      double * __restrict__ phi_z_g = phi_z[group];

      for(int d = 0; d < num_local_directions; d++){
        double **ell_plus_d = ell_plus[d];
        double psi_z_g_d = 0.0;

        for(int n = 0; n < num_moments; n++){
          double *ell_plus_d_n = ell_plus_d[n];
          double *phi_z_g_n = phi_z_g + n*(n+1)/2;

          for(int m = 0; m <= n; m++){
            psi_z_g_d += ell_plus_d_n[m] * phi_z_g_n[m];
          }
        }

        psi_z_g[d] = psi_z_g_d;
      }
    }
  }
}


static void zgd_nmd_LPlusTimes_3d(zgd_nmd_Param &p) {
  // Grab parameters
  double ***psi = p.psi.data;
  double ***phi = p.phi.data;
  double ***ell_plus = p.ell_plus.data;

  int num_groups = p.num_groups;
  int num_local_directions = p.num_directions;
  int num_zones = p.num_zones;

  int num_moments = p.num_moments;

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(int z = 0; z < num_zones; z++){
    double **psi_z = psi[z];
    double **phi_z = phi[z];

    for(int group = 0; group < num_groups; ++group){
      double * __restrict__ psi_z_g = psi_z[group];
      double * __restrict__ phi_z_g = phi_z[group];

      for(int d = 0; d < num_local_directions; d++){
        double **ell_plus_d = ell_plus[d];
        double psi_z_g_d = 0.0;

        for(int n = 0; n < num_moments; n++){
          double * __restrict__ ell_plus_d_n = ell_plus_d[n];
          double * __restrict__ phi_z_g_n = phi_z_g + n*n;

          int mmax  = 2*n+1;
          for(int m = 0; m < mmax; m++){
            psi_z_g_d += ell_plus_d_n[m] * phi_z_g_n[m];
          }
        }

        psi_z_g[d] = psi_z_g_d;
      }
    }
  }
}


/**
 * LPlusTimes routine with data nesting of:
 *    Psi[group][direction][zone]
 *    Phi[moment][zone]
 */
void zgd_nmd_Param::LPlusTimes(void) {

  if(geometry_type == 1){
    zgd_nmd_LPlusTimes_1d(*this);
  }
  else if(geometry_type == 2){
    zgd_nmd_LPlusTimes_2d(*this);
  }
  else {
    zgd_nmd_LPlusTimes_3d(*this);
  }
}
