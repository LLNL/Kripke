#include<Kripke/Kernel/Kernel_3d_GDZ.h>
#include<Kripke/User_Data.h>
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

void Kernel_3d_GDZ::LTimes(Grid_Data *grid_data) {
  // Outer parameters
  int num_zones = grid_data->num_zones;
  int nidx = grid_data->nm_table.size();
	int num_directions = grid_data->ell->directions;
	
  // Clear phi
  grid_data->phi->clear(0.0);

  // Loop over Group Sets
  int num_group_sets = grid_data->gd_sets.size();
  for (int gset = 0; gset < num_group_sets; ++gset) {
    std::vector<Group_Dir_Set> &dir_sets = grid_data->gd_sets[gset];
    int num_dir_sets = dir_sets.size();

    // Loop over Direction Sets
    for (int dset = 0; dset < num_dir_sets; ++dset) {
      Group_Dir_Set &gd_set = dir_sets[dset];

      // Get dimensioning
      int num_local_groups = gd_set.num_groups;
      int group0 = gd_set.group0;
      int num_local_directions = gd_set.num_directions;
      int dir0 = gd_set.direction0;

      /* 3D Cartesian Geometry */
      double * KRESTRICT phi = grid_data->phi->ptr(group0, 0, 0);
      double * KRESTRICT psi_ptr = gd_set.psi->ptr();

      for (int group = 0; group < num_local_groups; ++group) {
        double * KRESTRICT ell_nm = grid_data->ell->ptr(0, dir0, 0);

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
          ell_nm += num_directions;
          phi += num_zones;
        }
        psi_ptr += num_zones*num_local_directions;
      }
    } // Direction Set
  } // Group Set
}

void Kernel_3d_GDZ::LPlusTimes(Grid_Data *grid_data) {
  // Outer parameters
  int num_zones = grid_data->num_zones;
  int nidx = grid_data->nm_table.size();

  // Loop over Group Sets
  int num_group_sets = grid_data->gd_sets.size();
  for (int gset = 0; gset < num_group_sets; ++gset) {
    std::vector<Group_Dir_Set> &dir_sets = grid_data->gd_sets[gset];
    int num_dir_sets = dir_sets.size();

    // Loop over Direction Sets
    for (int dset = 0; dset < num_dir_sets; ++dset) {
      Group_Dir_Set &gd_set = dir_sets[dset];

      // Get dimensioning
      int num_local_groups = gd_set.num_groups;
      int group0 = gd_set.group0;
      int num_local_directions = gd_set.num_directions;
      int dir0 = gd_set.direction0;

      // Get Variables
      gd_set.rhs->clear(0.0);

      /* 3D Cartesian Geometry */
      double *ell_plus_ptr = grid_data->ell_plus->ptr(0, dir0, 0);

      double * KRESTRICT phi_out_ptr = grid_data->phi_out->ptr(group0, 0, 0);
      double * KRESTRICT rhs = gd_set.rhs->ptr();

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
    } // Direction Set
  } // Group Set
}

/* Sweep routine for Diamond-Difference */
/* Macros for offsets with fluxes on cell faces */
#define I_PLANE_INDEX(j, k) (k)*(local_jmax) + (j)
#define J_PLANE_INDEX(i, k) (k)*(local_imax) + (i)
#define K_PLANE_INDEX(i, j) (j)*(local_imax) + (i)
#define Zonal_INDEX(i, j, k) (i) + (local_imax)*(j) \
  + (local_imax)*(local_jmax)*(k)

void Kernel_3d_GDZ::sweep(Grid_Data *grid_data, Group_Dir_Set *gd_set,
    double *i_plane_ptr, double *j_plane_ptr, double *k_plane_ptr) {
  int num_directions = gd_set->num_directions;
  int num_groups = gd_set->num_groups;
  int num_zones = grid_data->num_zones;

  Directions *direction = gd_set->directions;

  int local_imax = grid_data->nzones[0];
  int local_jmax = grid_data->nzones[1];
  int local_kmax = grid_data->nzones[2];

  double * dx = &grid_data->deltas[0][0];
  double * dy = &grid_data->deltas[1][0];
  double * dz = &grid_data->deltas[2][0];

  // Alias the MPI data with a SubTVec for the face data
  SubTVec i_plane(nestingPsi(), num_groups, num_directions,
      local_jmax * local_kmax, i_plane_ptr);
  SubTVec j_plane(nestingPsi(), num_groups, num_directions,
      local_imax * local_kmax, j_plane_ptr);
  SubTVec k_plane(nestingPsi(), num_groups, num_directions,
      local_imax * local_jmax, k_plane_ptr);

  // All directions have same id,jd,kd, since these are all one Direction Set
  // So pull that information out now
  int octant = direction[0].octant;
  Grid_Sweep_Block const &extent = grid_data->octant_extent[octant];

#ifdef KRIPKE_USE_OPENMP
#pragma omp parallel for
#endif
  for (int group = 0; group < num_groups; ++group) {

    std::vector<double> xcos_dxi_all(local_imax);
    std::vector<double> ycos_dyj_all(local_jmax);
    std::vector<double> zcos_dzk_all(local_kmax);
    double * KRESTRICT sigt_g = grid_data->sigt->ptr(group+gd_set->group0, 0, 0);

    for (int d = 0; d < num_directions; ++d) {
      double * KRESTRICT psi_g_d = gd_set->psi->ptr(group, d, 0);
      double * KRESTRICT rhs_g_d = gd_set->rhs->ptr(group, d, 0);
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

