#include "../Param.h"

static void zdg_nmd_LPlusTimes_1d(zdg_nmd_Param &p) {
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

    for(int d = 0; d < num_local_directions; d++){
      double **ell_plus_d = ell_plus[d];
      double * __restrict__ psi_z_d = psi_z[d];

      for(int group = 0; group < num_groups; ++group){
        psi_z_d[group] = 0.0;
      }

      for(int n = 0; n < num_moments; n++){
        double *ell_plus_d_n = ell_plus_d[n];
        double ell_plus_d_n_m = ell_plus_d_n[0]; // m = 0;
        double * __restrict__ phi_z_n = phi_z[n];

        for(int group = 0; group < num_groups; ++group){
          psi_z_d[group] += ell_plus_d_n_m * phi_z_n[group];
        }
      }
    }
  }

}


static void zdg_nmd_LPlusTimes_2d(zdg_nmd_Param &p) {
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

    for(int d = 0; d < num_local_directions; d++){
      double **ell_plus_d = ell_plus[d];
      double * __restrict__ psi_z_d = psi_z[d];

      for(int group = 0; group < num_groups; ++group){
        psi_z_d[group] = 0.0;
      }

      for(int n = 0; n < num_moments; n++){
        double *ell_plus_d_n = ell_plus_d[n];

        double **phi_z_n = phi_z + n*(n+1)/2;

        for(int m = 0; m <= n; m++){
          double * __restrict__ phi_z_n_m = phi_z_n[m];
          double ell_plus_d_n_m = ell_plus_d_n[m];
          for(int group = 0; group < num_groups; ++group){
            psi_z_d[group] += ell_plus_d_n_m * phi_z_n_m[group];
          }
        }
      }
    }
  }
}


static void zdg_nmd_LPlusTimes_3d(zdg_nmd_Param &p) {
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

    for(int d = 0; d < num_local_directions; d++){
      double **ell_plus_d = ell_plus[d];
      double * __restrict__ psi_z_d = psi_z[d];

      for(int group = 0; group < num_groups; ++group){
        psi_z_d[group] = 0.0;
      }

      for(int n = 0; n < num_moments; n++){
        double *ell_plus_d_n = ell_plus_d[n];

        double **phi_z_n = phi_z + n*n;

        for(int m = -n; m <= n; m++){
          double * __restrict__ phi_z_n_m = phi_z_n[m+n];
          double ell_plus_d_n_m = ell_plus_d_n[m+n];
          for(int group = 0; group < num_groups; ++group){
            psi_z_d[group] += ell_plus_d_n_m * phi_z_n_m[group];
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
void zdg_nmd_Param::LPlusTimes(void) {

  if(geometry_type == 1){
    zdg_nmd_LPlusTimes_1d(*this);
  }
  else if(geometry_type == 2){
    zdg_nmd_LPlusTimes_2d(*this);
  }
  else {
    zdg_nmd_LPlusTimes_3d(*this);
  }
}
