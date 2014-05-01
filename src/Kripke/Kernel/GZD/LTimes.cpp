#include "../Param.h"

static void gzd_nmd_LTimes_1d(gzd_nmd_Param &p) {
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
    for(int z = 0; z < num_zones; z++){
      double * __restrict__ psi_g_z = psi_zonal[z];
      double * __restrict__ phi_g_z = phi_g[z];

      for(int n = 0; n < num_moments; n++){
        double phi_g_z_n = 0.0;
        double * __restrict__ ell_n_m = ell[n][0]; // m == 0

        for(int d = 0; d < num_local_directions; d++){
          double ell_n_m_d = ell_n_m[d];
          double psi_g_z_d = psi_g_z[d];
          phi_g_z_n += ell_n_m_d * psi_g_z_d;
        }
        phi_g_z[n] = phi_g_z_n;
      }
    }
  }
}


static void gzd_nmd_LTimes_2d(gzd_nmd_Param &p) {
  // Grab parameters
  double ***psi = p.psi.data;
  double ***phi = p.phi.data;
  double ***ell = p.ell.data;

  int num_groups = p.num_groups;
  int num_local_directions = p.num_directions;
  int num_zones = p.num_zones;

  int num_moments = p.num_moments;

  /* 2D rho-z Cylindrical Geometry */
  for(int group = 0; group < num_groups; ++group){
    double **psi_zonal = psi[group];
    double **phi_g = phi[group];

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(int z = 0; z < num_zones; z++){
      double * __restrict__ psi_g_z = psi_zonal[z];
      double * __restrict__ phi_g_z = phi_g[z];

      for(int n = 0; n < num_moments; n++){
        double * __restrict__ phi_g_z_n = phi_g_z + (n*(n+1)/2);
        double **ell_n = ell[n];
        for(int m = 0; m <= n; m++){

          double * __restrict__ ell_n_m = ell_n[m];
          double phi_g_z_nm = 0.0;

          for(int d = 0; d < num_local_directions; d++){
            double ell_n_m_d = ell_n_m[d];
            double psi_g_z_d = psi_g_z[d];
            phi_g_z_nm += ell_n_m_d * psi_g_z_d;
          }

          phi_g_z_n[m] = phi_g_z_nm;
        }
      }
    }
  }
}

#if 1
static void gzd_nmd_LTimes_3d(gzd_nmd_Param &p) {
  // Grab parameters
  double ***psi = p.psi.data;
  double ***phi = p.phi.data;
  double ***ell = p.ell.data;

  int num_groups = p.num_groups;
  int num_local_directions = p.num_directions;
  int num_zones = p.num_zones;

  int num_moments = p.num_moments;

  /* 3D Cartesian Geometry */
  for(int group = 0; group < num_groups; ++group){
    double **psi_zonal = psi[group];
    double **phi_g = phi[group];

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(int z = 0; z < num_zones; z++){
      double * __restrict__ psi_g_z = psi_zonal[z];
      double * __restrict__ phi_g_z = phi_g[z];

      for(int n = 0; n < num_moments; n++){
        double * __restrict__ phi_g_z_n = phi_g_z + (n*n);
        double **ell_n = ell[n];
        int nn = 2*n;
        for(int m = 0; m <= nn; m++){

          double* __restrict__ ell_n_m = ell_n[m];
          double phi_g_z_nm = 0.0;

          for(int d = 0; d < num_local_directions; d++){
            double ell_n_m_d = ell_n_m[d];
            double psi_g_z_d = psi_g_z[d];
            phi_g_z_nm += ell_n_m_d * psi_g_z_d;
          }

          phi_g_z_n[m] = phi_g_z_nm;
        }
      }
    }
  }
}

#else

static void gzd_nmd_LTimes_3d_inner(gzd_nmd_Param &p) {
  // Grab parameters
  double ***psi = p.psi.data;
  double ***phi = p.phi.data;
  double ***ell = p.ell.data;

  int num_groups = p.num_groups;
  int num_local_directions = p.num_directions;
  int num_zones = p.num_zones;

  int num_moments = p.num_moments;

  /* 3D Cartesian Geometry */
  for(int group = 0; group < num_groups; ++group){
    double **psi_zonal = psi[group];
    double **phi_g = phi[group];

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(int z = 0; z < num_zones; z++){
      double * __restrict__ psi_g_z = psi_zonal[z];
      double * __restrict__ phi_g_z = phi_g[z];

      for(int n = 0; n < num_moments; n++){
        double * __restrict__ phi_g_z_n = phi_g_z + (n*n);
        double **ell_n = ell[n];
        int nn = 2*n;
        for(int m = 0; m <= nn; m++){

          double* __restrict__ ell_n_m = ell_n[m];
          double phi_g_z_nm = 0.0;

          for(int d = 0; d < num_local_directions; d++){
            double ell_n_m_d = ell_n_m[d];
            double psi_g_z_d = psi_g_z[d];
            phi_g_z_nm += ell_n_m_d * psi_g_z_d;
          }

          phi_g_z_n[m] = phi_g_z_nm;
        }
      }
    }
  }
}

static void gzd_nmd_LTimes_3d_outer(gzd_nmd_Param &p, int G, int D, int Z) {
  int num_groups = p.num_groups;
  int num_local_directions = p.num_directions;
  int num_zones = p.num_zones;

  D >
}

static void gzd_nmd_LTimes_3d(gzd_nmd_Param &p) {
  int num_groups = p.num_groups;
  int num_local_directions = p.num_directions;
  int num_zones = p.num_zones;

  gzd_nmd_LTimes_3d_outer(p, num_groups, num_local_directions, num_zones);
}


#endif

/**
 * LTimes routine with data nesting of:
 *    Psi[group][direction][zone]
 *    Phi[moment][zone]
 */
void gzd_nmd_Param::LTimes(void) {

  if(geometry_type == 1){
    gzd_nmd_LTimes_1d(*this);
  }
  else if(geometry_type == 2){
    gzd_nmd_LTimes_2d(*this);
  }
  else {
    gzd_nmd_LTimes_3d(*this);
  }
}
