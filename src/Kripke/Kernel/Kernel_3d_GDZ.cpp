#include<Kripke/Kernel/Kernel_3d_GDZ.h>
#include<Kripke/Grid.h>
#include<Kripke/SubTVec.h>

Kernel_3d_GDZ::Kernel_3d_GDZ() {

}

Kernel_3d_GDZ::~Kernel_3d_GDZ() {

}

Nesting_Order Kernel_3d_GDZ::nestingPsi(void) const {
  return NEST_GDZ;
}

Nesting_Order Kernel_3d_GDZ::nestingPhi(void) const {
  return NEST_GDZ;
}

Nesting_Order Kernel_3d_GDZ::nestingSigt(void) const {
  return NEST_DGZ;
}

Nesting_Order Kernel_3d_GDZ::nestingEll(void) const {
  return NEST_ZGD;
}

Nesting_Order Kernel_3d_GDZ::nestingEllPlus(void) const {
  return NEST_ZDG;
}

void Kernel_3d_GDZ::LTimes(Grid_Data *grid_data) {
  // Outer parameters
  int nidx = grid_data->total_num_moments;

  // Clear phi
  for(int ds = 0;ds < grid_data->num_zone_sets;++ ds){
    grid_data->phi[ds]->clear(0.0);
  }

  // Loop over Subdomains
  int num_subdomains = grid_data->subdomains.size();
  for (int sdom_id = 0; sdom_id < num_subdomains; ++ sdom_id){
    Subdomain &sdom = grid_data->subdomains[sdom_id];

    // Get dimensioning
    int num_zones = sdom.num_zones;
    int num_local_groups = sdom.num_groups;
    int group0 = sdom.group0;
    int num_local_directions = sdom.num_directions;

    /* 3D Cartesian Geometry */
    double * KRESTRICT phi = sdom.phi->ptr(group0, 0, 0);
    double * KRESTRICT psi_ptr = sdom.psi->ptr();

    for (int group = 0; group < num_local_groups; ++group) {
      double * KRESTRICT ell_nm = sdom.ell->ptr();

      for(int nm_offset = 0;nm_offset < nidx;++nm_offset){
        double * KRESTRICT psi = psi_ptr;

        for (int d = 0; d < num_local_directions; d++) {
          double ell_nm_d = ell_nm[d];

#ifdef KRIPKE_USE_OPENMP
#pragma omp parallel for
#endif
          for(int z = 0;z < num_zones; ++ z){
            phi[z] += ell_nm_d * psi[z];
          }

          psi += num_zones;
        }
        ell_nm += num_local_directions;
        phi += num_zones;
      }
      psi_ptr += num_zones*num_local_directions;
    }
  } // Subdomain
}

void Kernel_3d_GDZ::LPlusTimes(Grid_Data *grid_data) {
  // Outer parameters
  int nidx = grid_data->total_num_moments;

  // Loop over Subdomains
  int num_subdomains = grid_data->subdomains.size();
  for (int sdom_id = 0; sdom_id < num_subdomains; ++ sdom_id){
    Subdomain &sdom = grid_data->subdomains[sdom_id];

    // Get dimensioning
    int num_zones = sdom.num_zones;
    int num_local_groups = sdom.num_groups;
    int group0 = sdom.group0;
    int num_local_directions = sdom.num_directions;

    // Get Variables
    sdom.rhs->clear(0.0);

    /* 3D Cartesian Geometry */
    double *ell_plus_ptr = sdom.ell_plus->ptr();

    double * KRESTRICT phi_out_ptr = sdom.phi_out->ptr(group0, 0, 0);
    double * KRESTRICT rhs = sdom.rhs->ptr();

    for (int group = 0; group < num_local_groups; ++group) {
      double *ell_plus = ell_plus_ptr;

      for (int d = 0; d < num_local_directions; d++) {
        double * KRESTRICT phi_out = phi_out_ptr;

        for(int nm_offset = 0;nm_offset < nidx;++nm_offset){
          double ell_plus_d_nm = ell_plus[nm_offset];

#ifdef KRIPKE_USE_OPENMP
#pragma omp parallel for
#endif
          for(int z = 0;z < num_zones; ++ z){
            rhs[z] += ell_plus_d_nm * phi_out[z];
          }
          phi_out += num_zones;
        }
        ell_plus += nidx;
        rhs += num_zones;
      }
      phi_out_ptr += num_zones*nidx;
    }
  } // Subdomain
}

/* Sweep routine for Diamond-Difference */
/* Macros for offsets with fluxes on cell faces */
#define I_PLANE_INDEX(j, k) (k)*(local_jmax) + (j)
#define J_PLANE_INDEX(i, k) (k)*(local_imax) + (i)
#define K_PLANE_INDEX(i, j) (j)*(local_imax) + (i)
#define Zonal_INDEX(i, j, k) (i) + (local_imax)*(j) \
  + (local_imax)*(local_jmax)*(k)

void Kernel_3d_GDZ::sweep(Subdomain *sdom) {
  int num_directions = sdom->num_directions;
  int num_groups = sdom->num_groups;
  int num_zones = sdom->num_zones;

  Directions *direction = sdom->directions;

  int local_imax = sdom->nzones[0];
  int local_jmax = sdom->nzones[1];
  int local_kmax = sdom->nzones[2];

  double * dx = &sdom->deltas[0][0];
  double * dy = &sdom->deltas[1][0];
  double * dz = &sdom->deltas[2][0];

  // Upwind/Downwind face flux data
  SubTVec &i_plane = *sdom->plane_data[0];
  SubTVec &j_plane = *sdom->plane_data[1];
  SubTVec &k_plane = *sdom->plane_data[2];

  // All directions have same id,jd,kd, since these are all one Direction Set
  // So pull that information out now
  Grid_Sweep_Block const &extent = sdom->sweep_block;

#ifdef KRIPKE_USE_OPENMP
#pragma omp parallel for
#endif
  for (int group = 0; group < num_groups; ++group) {

    std::vector<double> xcos_dxi_all(local_imax);
    std::vector<double> ycos_dyj_all(local_jmax);
    std::vector<double> zcos_dzk_all(local_kmax);
    double * KRESTRICT sigt_g = sdom->sigt->ptr(group, 0, 0);

    for (int d = 0; d < num_directions; ++d) {
      double * KRESTRICT psi_g_d = sdom->psi->ptr(group, d, 0);
      double * KRESTRICT rhs_g_d = sdom->rhs->ptr(group, d, 0);
      double * KRESTRICT i_plane_g_d = i_plane.ptr(group, d, 0);
      double * KRESTRICT j_plane_g_d = j_plane.ptr(group, d, 0);
      double * KRESTRICT k_plane_g_d = k_plane.ptr(group, d, 0);

      double xcos = direction[d].xcos;
      double ycos = direction[d].ycos;
      double zcos = direction[d].zcos;
      for (int i = 0; i < local_imax; ++i) {
        double dxi = dx[i + 1];
        xcos_dxi_all[i] = 2.0 * xcos / dxi;
      }

      for (int j = 0; j < local_jmax; ++j) {
        double dyj = dy[j + 1];
        ycos_dyj_all[j] = 2.0 * ycos / dyj;
      }

      for (int k = 0; k < local_kmax; ++k) {
        double dzk = dz[k + 1];
        zcos_dzk_all[k] = 2.0 * zcos / dzk;
      }

      /*  Perform transport sweep of the grid 1 cell at a time.   */
      for (int k = extent.start_k; k != extent.end_k; k += extent.inc_k) {
        double zcos_dzk = zcos_dzk_all[k];
        for (int j = extent.start_j; j != extent.end_j; j += extent.inc_j) {
          double ycos_dyj = ycos_dyj_all[j];
          int z_idx = Zonal_INDEX(extent.start_i, j, k);
          for (int i = extent.start_i; i != extent.end_i; i += extent.inc_i) {
            double xcos_dxi = xcos_dxi_all[i];

            /* Calculate new zonal flux */
            double psi_g_d_z = (rhs_g_d[z_idx]
                + i_plane_g_d[I_PLANE_INDEX(j, k)] * xcos_dxi
                + j_plane_g_d[J_PLANE_INDEX(i, k)] * ycos_dyj
                + k_plane_g_d[K_PLANE_INDEX(i, j)] * zcos_dzk)
                / (xcos_dxi + ycos_dyj + zcos_dzk
                    + sigt_g[z_idx]);
            psi_g_d[z_idx] = psi_g_d_z;

            /* Apply diamond-difference relationships */
            i_plane_g_d[I_PLANE_INDEX(j, k)] = 2.0 * psi_g_d_z
                - i_plane_g_d[I_PLANE_INDEX(j, k)];
            j_plane_g_d[J_PLANE_INDEX(i, k)] = 2.0 * psi_g_d_z
                - j_plane_g_d[J_PLANE_INDEX(i, k)];
            k_plane_g_d[K_PLANE_INDEX(i, j)] = 2.0 * psi_g_d_z
                - k_plane_g_d[K_PLANE_INDEX(i, j)];

            z_idx += extent.inc_i;
          }
        }
      }
    }
  }

}

