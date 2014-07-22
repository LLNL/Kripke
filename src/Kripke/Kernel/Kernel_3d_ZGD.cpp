#include<Kripke/Kernel/Kernel_3d_ZGD.h>
#include<Kripke/User_Data.h>
#include<Kripke/SubTVec.h>

Kernel_3d_ZGD::Kernel_3d_ZGD() {

}

Kernel_3d_ZGD::~Kernel_3d_ZGD() {

}

Nesting_Order Kernel_3d_ZGD::nestingPsi(void) const {
  return NEST_ZGD;
}

Nesting_Order Kernel_3d_ZGD::nestingPhi(void) const {
  return NEST_ZGD;
}

void Kernel_3d_ZGD::LTimes(Grid_Data *grid_data) {
  // Outer parameters
  int num_zones = grid_data->num_zones;
  int num_moments = grid_data->num_moments;
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
      double * KRESTRICT ell_ptr = grid_data->ell->ptr(0, dir0, 0);

#ifdef KRIPKE_USE_OPENMP
#pragma omp parallel for
#endif
      for(int z = 0;z < num_zones;++ z){
        double * KRESTRICT psi = gd_set.psi->ptr(0, 0, z);
        double * KRESTRICT phi = grid_data->phi->ptr(group0, 0, z);
        for(int group = 0;group < num_local_groups;++ group){
          double * KRESTRICT ell_d = ell_ptr;

          for (int d = 0; d < num_local_directions; d++) {
            double psi_d = psi[d];

            for(int nm_offset = 0;nm_offset < nidx;++nm_offset){
              phi[nm_offset] += ell_d[nm_offset] * psi_d;
            }
            ell_d += nidx;
          }

          psi += num_local_directions;
          phi += nidx;
        }
      }

    } // Direction Set
  } // Group Set
}

void Kernel_3d_ZGD::LPlusTimes(Grid_Data *grid_data) {
  // Outer parameters
  int num_zones = grid_data->num_zones;
  int num_moments = grid_data->num_moments;
  int nidx = grid_data->nm_table.size();
  int num_directions = grid_data->ell_plus->directions;

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
      double * KRESTRICT ell_plus_ptr = grid_data->ell_plus->ptr(0, dir0, 0);


#ifdef KRIPKE_USE_OPENMP
#pragma omp parallel for
#endif
      for(int z = 0;z < num_zones;++ z){
        double * KRESTRICT rhs = gd_set.rhs->ptr(0, 0, z);
        double * KRESTRICT phi_out = grid_data->phi_out->ptr(group0,0, z);
        for(int group = 0;group < num_local_groups;++ group){
          double * KRESTRICT ell_plus_d = ell_plus_ptr;

          for (int d = 0; d < num_local_directions; d++) {
            double rhs_acc = 0.0;

            for(int nm_offset = 0;nm_offset < nidx;++nm_offset){
               rhs_acc += ell_plus_d[nm_offset] * phi_out[nm_offset];
            }
            rhs[d] += rhs_acc;

            ell_plus_d += nidx;
          }
          rhs += num_local_directions;
          phi_out += nidx;
        }
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

void Kernel_3d_ZGD::sweep(Grid_Data *grid_data, Group_Dir_Set *gd_set,
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

  /*  Perform transport sweep of the grid 1 cell at a time.   */

  for (int k = extent.start_k; k != extent.end_k; k += extent.inc_k) {
    double dzk = dz[k + 1];
    double two_dz = 2.0 / dzk;
    for (int j = extent.start_j; j != extent.end_j; j += extent.inc_j) {
      double dyj = dy[j + 1];
      double two_dy = 2.0 / dyj;
      for (int i = extent.start_i; i != extent.end_i; i += extent.inc_i) {
        double dxi = dx[i + 1];
        double two_dx = 2.0 / dxi;

        int z = Zonal_INDEX(i, j, k);
        double * sigt_z = grid_data->sigt->ptr(gd_set->group0, 0, z);
#ifdef KRIPKE_USE_OPENMP
#pragma omp parallel for
#endif
        for (int group = 0; group < num_groups; ++group) {
          double * KRESTRICT psi_z_g = gd_set->psi->ptr(group, 0, z);
          double * KRESTRICT rhs_z_g = gd_set->rhs->ptr(group, 0, z);

          double * KRESTRICT psi_lf_z_g = i_plane.ptr(group, 0, I_PLANE_INDEX(j, k));
          double * KRESTRICT psi_fr_z_g = j_plane.ptr(group, 0, J_PLANE_INDEX(i, k));
          double * KRESTRICT psi_bo_z_g = k_plane.ptr(group, 0, K_PLANE_INDEX(i, j));

          for (int d = 0; d < num_directions; ++d) {

            double xcos = direction[d].xcos;
            double ycos = direction[d].ycos;
            double zcos = direction[d].zcos;

            double zcos_dzk = zcos * two_dz;
            double ycos_dyj = ycos * two_dy;
            double xcos_dxi = xcos * two_dx;

            /* Calculate new zonal flux */
            double psi_z_g_d = (rhs_z_g[d]
                + psi_lf_z_g[d] * xcos_dxi
                + psi_fr_z_g[d] * ycos_dyj
                + psi_bo_z_g[d] * zcos_dzk)
                / (xcos_dxi + ycos_dyj + zcos_dzk + sigt_z[group]);

            psi_z_g[d] = psi_z_g_d;

            /* Apply diamond-difference relationships */
            psi_lf_z_g[d] = 2.0 * psi_z_g_d - psi_lf_z_g[d];
            psi_fr_z_g[d] = 2.0 * psi_z_g_d - psi_fr_z_g[d];
            psi_bo_z_g[d] = 2.0 * psi_z_g_d - psi_bo_z_g[d];
          }
        }
      }
    }
  }
}

