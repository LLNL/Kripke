#include "../Param.h"

static void zdg_nmd_LTimes_1d(zdg_nmd_Param &p) {
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

    for(int n = 0; n < num_moments; n++){
      double * __restrict__ ell_n_m = ell[n][0];   // m == 0
      double * __restrict__ phi_z_n = phi_z[n];

      for(int group = 0; group < num_groups; ++group){
        phi_z_n[group] = 0;
      }

      for(int d = 0; d < num_local_directions; d++){
        double * __restrict__ psi_z_d = psi_z[d];
        double ell_n_m_d = ell_n_m[d];

        for(int group = 0; group < num_groups; ++group){
          double psi_z_d_g = psi_z_d[group];
          phi_z_n[group] += ell_n_m_d * psi_z_d_g;
        }
      }
    }
  }

}


static void zdg_nmd_LTimes_2d(zdg_nmd_Param &p) {
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

    for(int n = 0; n < num_moments; n++){
      double ** __restrict__ phi_z_n = phi_z +(n*(n+1)/2);
      double **ell_n = ell[n];
      for(int m = 0; m <= n; m++){
        double * __restrict__ ell_n_m = ell_n[m];
        double * __restrict__ phi_z_nm = phi_z_n[m];

        for(int group = 0; group < num_groups; ++group){
          phi_z_nm[group] = 0;
        }

        for(int d = 0; d < num_local_directions; d++){
          double * __restrict__ psi_z_d = psi_z[d];
          double ell_n_m_d = ell_n_m[d];

          for(int group = 0; group < num_groups; ++group){
            double psi_z_d_g = psi_z_d[group];
            phi_z_nm[group] += ell_n_m_d * psi_z_d_g;
          }
        }
      }
    }
  }
}


static void zdg_nmd_LTimes_3d(zdg_nmd_Param &p) {
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

    for(int n = 0; n < num_moments; n++){
      double ** __restrict__ phi_z_n = phi_z + n*n + n;
      double **ell_n = ell[n];
      for(int m = -n; m <= n; m++){
        double * __restrict__ ell_n_m = ell_n[m+n];
        double * __restrict__ phi_z_nm = phi_z_n[m];

        for(int group = 0; group < num_groups; ++group){
          phi_z_nm[group] = 0;
        }

        for(int d = 0; d < num_local_directions; d++){
          double * __restrict__ psi_z_d = psi_z[d];
          double ell_n_m_d = ell_n_m[d];

          for(int group = 0; group < num_groups; ++group){
            double psi_z_d_g = psi_z_d[group];
            phi_z_nm[group] += ell_n_m_d * psi_z_d_g;
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
void zdg_nmd_Param::LTimes(void) {

  if(geometry_type == 1){
    zdg_nmd_LTimes_1d(*this);
  }
  else if(geometry_type == 2){
    zdg_nmd_LTimes_2d(*this);
  }
  else {
    zdg_nmd_LTimes_3d(*this);
  }
}
