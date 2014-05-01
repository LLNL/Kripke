#include "../Param.h"

void gzd_nmd_LPlusTimes_1d(gzd_nmd_Param &p) {
  // Grab parameters
  double ***psi = p.psi.data;
  double ***phi = p.phi.data;
  double ***ell_plus = p.ell_plus.data;

  int num_groups = p.num_groups;
  int num_local_directions = p.num_directions;
  int num_zones = p.num_zones;

  int num_moments = p.num_moments;

  /* 1D Spherical or Slab Geometry */

  for(int group = 0; group < num_groups; ++group){
    double **phi_g = phi[group];
    double **psi_g = psi[group];

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(int i=0; i < num_zones; i++){
      double const * __restrict__ phi_g_z = phi_g[i];
      double * __restrict__ psi_g_z = psi_g[i];

      for(int d=0; d < num_local_directions; d++){
        //double **ell_plus_d = ell_plus[d];
        double const * __restrict__ ell_plus_d_m = ell_plus[d][0];
        double psi_g_z_d = 0.0;
        for(int n=0; n< num_moments; n++){
          //double *ell_plus_d_n = ell_plus_d[n];
          //double ell_plus_d_n_m = ell_plus_d_n[0]; // m = 0;
          //psi_g_z_d += ell_plus_d_n_m * phi_g_z[n];
          psi_g_z_d += ell_plus_d_m[n] * phi_g_z[n];
        }
        psi_g_z[d] = psi_g_z_d;
      }
    }
  }
}


void gzd_nmd_LPlusTimes_2d(gzd_nmd_Param &p) {
  // Grab parameters
  double ***psi = p.psi.data;
  double ***phi = p.phi.data;
  double ***ell_plus = p.ell_plus.data;

  int num_groups = p.num_groups;
  int num_local_directions = p.num_directions;
  int num_zones = p.num_zones;

  int num_moments = p.num_moments;

  /* 2D rho-z Cylindrical Geometry */

  for(int group = 0; group < num_groups; ++group){
    double **phi_g = phi[group];
    double **psi_g = psi[group];

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(int i=0; i < num_zones; i++){
      double * __restrict__ phi_g_z = phi_g[i];
      double * __restrict__ psi_g_z = psi_g[i];
      for(int d=0; d < num_local_directions; d++){
        double **ell_plus_d = ell_plus[d];
        double psi_g_z_d = 0.0;

        for(int n=0; n< num_moments; n++){
          double * __restrict__ ell_plus_d_n = ell_plus_d[n];
          double * __restrict__ phi_g_z_n = phi_g_z + n*(n+1)/2;

          for(int m = 0; m <= n; m++){
            psi_g_z_d += ell_plus_d_n[m] * phi_g_z_n[m];
          }
        }

        psi_g_z[d] = psi_g_z_d;
      }
    }
  }
}


void gzd_nmd_LPlusTimes_3d(gzd_nmd_Param &p) {
  // Grab parameters
  double ***psi = p.psi.data;
  double ***phi = p.phi.data;
  double ***ell_plus = p.ell_plus.data;

  int num_groups = p.num_groups;
  int num_local_directions = p.num_directions;
  int num_zones = p.num_zones;

  int num_moments = p.num_moments;

  int nmom = num_moments;
  nmom = (nmom-1)*(nmom-1)+2*(nmom-1)+1;

  /* 3D Cartesian Geometry */
  for(int group = 0; group < num_groups; ++group){
    double **phi_g = phi[group];
    double **psi_g = psi[group];

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(int i=0; i < num_zones; i++){
      double const * __restrict__ phi_g_z = phi_g[i];
      double * __restrict__ psi_g_z = psi_g[i];
      for(int d=0; d < num_local_directions; d++){
        double **ell_plus_d = ell_plus[d];
        double psi_g_z_d = 0.0;
#if 0
        for(int n=0; n< num_moments; n++){
          int nn = n*n;
          int n2 = 2*n;
          double const * __restrict__ ell_plus_d_n = ell_plus_d[n];
          double const * __restrict__ phi_g_z_n = phi_g_z + nn;

          double psi_g_z_d_m = 0.0;
          for(int m = 0; m <= n2; m++){
            psi_g_z_d_m += ell_plus_d_n[m] * phi_g_z_n[m];
          }
          psi_g_z_d += psi_g_z_d_m;
        }
#else

        double const * __restrict__ ell_plus_d_n = ell_plus_d[0];
        for(int nm = 0; nm < nmom; ++nm){
          psi_g_z_d += ell_plus_d_n[nm] * phi_g_z[nm];
        }
#endif
        psi_g_z[d] = psi_g_z_d;
      }
    }
  }
}


/**
 * LPlusTimes routine with data nesting of:
 *    Psi[group][direction][zone]
 *    Phi[moment][zone]
 */
void gzd_nmd_Param::LPlusTimes(void) {

  if(geometry_type == 1){
    gzd_nmd_LPlusTimes_1d(*this);
  }
  else if(geometry_type == 2){
    gzd_nmd_LPlusTimes_2d(*this);
  }
  else {
    gzd_nmd_LPlusTimes_3d(*this);
  }
}
