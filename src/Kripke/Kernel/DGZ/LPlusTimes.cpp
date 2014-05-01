#include "../Param.h"

static void dgz_nmd_LPlusTimes_1d(dgz_nmd_Param &p) {
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
      for(int group = 0; group < num_groups; ++group){
        double * __restrict__ psi_d_g = psi_d[group];
        double * __restrict__ phi_n_g = phi_n[group];

        for(int z = 0; z < num_zones; z++){
          psi_d_g[z] += ell_plus_d_n_m * phi_n_g[z];
        }
      }
    }
  }

}


static void dgz_nmd_LPlusTimes_2d(dgz_nmd_Param &p) {
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
        for(int group = 0; group < num_groups; ++group){
          double * __restrict__ psi_d_g = psi_d[group];
          double * __restrict__ phi_n_g = phi_n_m[group];

          for(int z = 0; z < num_zones; z++){
            psi_d_g[z] += ell_plus_d_n_m * phi_n_g[z];
          }
        }
      }
    }
  }
}


static void dgz_nmd_LPlusTimes_3d(dgz_nmd_Param &p) {
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
        for(int group = 0; group < num_groups; ++group){
          double * __restrict__ psi_d_g = psi_d[group];
          double * __restrict__ phi_n_m_g = phi_n_m[group];

          for(int z = 0; z < num_zones; z++){
            psi_d_g[z] += ell_plus_d_n_m * phi_n_m_g[z];
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
void dgz_nmd_Param::LPlusTimes(void) {

  if(geometry_type == 1){
    dgz_nmd_LPlusTimes_1d(*this);
  }
  else if(geometry_type == 2){
    dgz_nmd_LPlusTimes_2d(*this);
  }
  else {
    dgz_nmd_LPlusTimes_3d(*this);
  }
}
