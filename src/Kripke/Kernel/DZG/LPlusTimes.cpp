#include "../Param.h"

static void dzg_nmd_LPlusTimes_1d(dzg_nmd_Param &p) {
  // Grab parameters
  double ***psi = p.psi.data;
  double ***phi = p.phi.data;
  double ***ell_plus = p.ell_plus.data;

  int num_groups = p.num_groups;
  int num_local_directions = p.num_directions;
  int num_zones = p.num_zones;

  int num_moments = p.num_moments;

  p.psi.clear(0.0);
  for(int d = 0; d < num_local_directions; d++){
    double **psi_d = psi[d];
    double **ell_plus_d = ell_plus[d];

    for(int n = 0; n < num_moments; n++){
      double **phi_n = phi[n];
      double *ell_plus_d_n = ell_plus_d[n];
      double ell_plus_d_n_m = ell_plus_d_n[0]; // m = 0;

#ifdef _OPENMP
#pragma omp parallel for
#endif
      for(int z = 0; z < num_zones; z++){
        double * __restrict__ psi_d_z = psi_d[z];
        double * __restrict__ phi_n_z = phi_n[z];

        for(int group = 0; group < num_groups; ++group){
          psi_d_z[group] += ell_plus_d_n_m * phi_n_z[group];
        }
      }
    }
  }

}


static void dzg_nmd_LPlusTimes_2d(dzg_nmd_Param &p) {
  // Grab parameters
  double ***psi = p.psi.data;
  double ***phi = p.phi.data;
  double ***ell_plus = p.ell_plus.data;

  int num_groups = p.num_groups;
  int num_local_directions = p.num_directions;
  int num_zones = p.num_zones;

  int num_moments = p.num_moments;

  p.psi.clear(0.0);
  for(int d = 0; d < num_local_directions; d++){
    double **psi_d = psi[d];
    double **ell_plus_d = ell_plus[d];

    for(int n = 0; n < num_moments; n++){
      double ***phi_n = phi + n*(n+1)/2;
      double *ell_plus_d_n = ell_plus_d[n];

      for(int m = 0; m <= n; m++){
        double **phi_n_m = phi_n[m];
        double ell_plus_d_n_m = ell_plus_d_n[m];

#ifdef _OPENMP
#pragma omp parallel for
#endif
        for(int z = 0; z < num_zones; z++){
          double * __restrict__ psi_d_z = psi_d[z];
          double * __restrict__ phi_n_z = phi_n_m[z];

          for(int group = 0; group < num_groups; ++group){
            psi_d_z[group] += ell_plus_d_n_m * phi_n_z[group];
          }
        }
      }
    }
  }
}


static void dzg_nmd_LPlusTimes_3d(dzg_nmd_Param &p) {
  // Grab parameters
  double ***psi = p.psi.data;
  double ***phi = p.phi.data;
  double ***ell_plus = p.ell_plus.data;

  int num_groups = p.num_groups;
  int num_local_directions = p.num_directions;
  int num_zones = p.num_zones;

  int num_moments = p.num_moments;

  p.psi.clear(0.0);
  for(int d = 0; d < num_local_directions; d++){
    double **psi_d = psi[d];
    double **ell_plus_d = ell_plus[d];

    for(int n = 0; n < num_moments; n++){
      double ***phi_n = phi + n*n;
      double *ell_plus_d_n = ell_plus_d[n];

      for(int m = 0; m <= 2*n; m++){
        double **phi_n_m = phi_n[m];
        double ell_plus_d_n_m = ell_plus_d_n[m];

#ifdef _OPENMP
#pragma omp parallel for
#endif
        for(int z = 0; z < num_zones; z++){
          double * __restrict__ psi_d_z = psi_d[z];
          double * __restrict__ phi_n_m_z = phi_n_m[z];

          for(int group = 0; group < num_groups; ++group){
            psi_d_z[group] += ell_plus_d_n_m * phi_n_m_z[group];
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
void dzg_nmd_Param::LPlusTimes(void) {

  if(geometry_type == 1){
    dzg_nmd_LPlusTimes_1d(*this);
  }
  else if(geometry_type == 2){
    dzg_nmd_LPlusTimes_2d(*this);
  }
  else {
    dzg_nmd_LPlusTimes_3d(*this);
  }
}
