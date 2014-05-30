#include<Kripke/Kernel/Kernel_3d_DGZ.h>
#include<Kripke/User_Data.h>
#include<Kripke/SubTVec.h>
#include<Kripke/LMat.h>

Kernel_3d_DGZ::Kernel_3d_DGZ() {

}

Kernel_3d_DGZ::~Kernel_3d_DGZ() {

}

Nesting_Order Kernel_3d_DGZ::nestingPsi(void) const {
  return NEST_DGZ;
}

Nesting_Order Kernel_3d_DGZ::nestingPhi(void) const {
  return NEST_DGZ;
}


void Kernel_3d_DGZ::LTimes(Grid_Data *grid_data) {
  // Outer parameters
  double ***ell = grid_data->ell->data;
  int num_zones = grid_data->num_zones;
  int num_moments = grid_data->num_moments;
  int nidx = grid_data->nm_table.size();

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
      int num_groups_zones = num_local_groups*num_zones;

      /* 3D Cartesian Geometry */
#ifdef KRIPKE_USE_OPENMP
#pragma omp parallel for
#endif
      for(int idx = 0;idx < nidx;++idx){
        int n = grid_data->nm_table[idx];
        int m = idx - n*n - n;

        double *ell_n_m = ell[n][m + n];
        int nm_offset = n*n + n + m;

        double *psi = gd_set.psi->ptr();
        for (int d = 0; d < num_local_directions; d++) {

          double ell_n_m_d = ell_n_m[d + dir0];
          double * KRESTRICT phi = grid_data->phi->ptr(group0, nm_offset, 0);
          double * KRESTRICT psi_ptr = psi;

          for(int gz = 0;gz < num_groups_zones;++ gz){
            phi[gz] += ell_n_m_d * psi_ptr[gz];

          }
          psi += num_groups_zones;
        }
      }

    } // Direction Set
  } // Group Set
}

void Kernel_3d_DGZ::LPlusTimes(Grid_Data *grid_data) {
  // Outer parameters
  double ***ell_plus = grid_data->ell_plus->data;
  int num_zones = grid_data->num_zones;
  int num_moments = grid_data->num_moments;
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
      int num_groups_zones = num_local_groups*num_zones;


      /* 3D Cartesian Geometry */
#if 0

#ifdef KRIPKE_USE_OPENMP
#pragma omp parallel for
#endif
      for (int d = 0; d < num_local_directions; d++) {
        double **ell_plus_d = ell_plus[d + dir0];
        double * KRESTRICT rhs = gd_set.rhs->ptr(0, d, 0);

        for(int gz = 0;gz < num_groups_zones;++ gz){
          rhs[gz] = 0;
        }

        for (int n = 0; n < num_moments; n++) {
          double *ell_plus_d_n = ell_plus_d[n];

          for (int m = -n; m <= n; m++) {
            int nm_offset = n*n + n + m;
            double ell_plus_d_n_m = ell_plus_d_n[n+m];

            double * KRESTRICT phi_out = grid_data->phi_out->ptr(group0, nm_offset, 0);

            for(int gz = 0;gz < num_groups_zones;++ gz){
              rhs[gz] += ell_plus_d_n_m * phi_out[gz];
            }
          }
        }
      }
#else

      int blk_size = 512;
#ifdef KRIPKE_USE_OPENMP
#pragma omp parallel for
#endif
      for(int gz_start = 0;gz_start < num_groups_zones;gz_start += blk_size){
        for (int d = 0; d < num_local_directions; d++) {
          double **ell_plus_d = ell_plus[d + dir0];
          double * KRESTRICT rhs = gd_set.rhs->ptr(0, d, 0);


          int gz_end = std::min(gz_start + blk_size, num_groups_zones);

          for(int gz = gz_start;gz < gz_end;++ gz){
            rhs[gz] = 0;
          }

          for(int nm_offset = 0;nm_offset < nidx;++nm_offset){
            int n = grid_data->nm_table[nm_offset];
            int m = nm_offset - n*n - n;
            double ell_plus_d_n_m = ell_plus_d[n][n+m];
            double * KRESTRICT phi_out = grid_data->phi_out->ptr(group0, nm_offset, 0);

            for(int gz = gz_start;gz < gz_end;++ gz){
              rhs[gz] += ell_plus_d_n_m * phi_out[gz];
            }
          }

        }
      }

#endif




    } // Direction Set
  } // Group Set
}

/* Sweep routine for Diamond-Difference */
/* Macros for offsets with fluxes on cell faces */
#define Left_INDEX(i, j, k) (i) + (local_imax_1)*(j) \
  + (local_imax_1)*(local_jmax)*(k)
#define Front_INDEX(i, j, k) (i) + (local_imax)*(j) \
  + (local_imax)*(local_jmax_1)*(k)
#define Bottom_INDEX(i, j, k) (i) + (local_imax)*(j) \
  + (local_imax)*(local_jmax)*(k)
#define I_PLANE_INDEX(j, k) (k)*(local_jmax) + (j)
#define J_PLANE_INDEX(i, k) (k)*(local_imax) + (i)
#define K_PLANE_INDEX(i, j) (j)*(local_imax) + (i)

#define Ghost_INDEX(i, j, k) (i) + (local_imax_2)*(j) \
  + (local_imax_2)*(local_jmax_2)*(k)
#define Zonal_INDEX(i, j, k) (i) + (local_imax)*(j) \
  + (local_imax)*(local_jmax)*(k)

void Kernel_3d_DGZ::sweep(Grid_Data *grid_data, Group_Dir_Set *gd_set,
    double *i_plane_ptr, double *j_plane_ptr, double *k_plane_ptr) {
  int num_directions = gd_set->num_directions;
  int num_groups = gd_set->num_groups;
  int num_zones = grid_data->num_zones;

  Directions *direction = gd_set->directions;

  int local_imax = grid_data->nzones[0];
  int local_jmax = grid_data->nzones[1];
  int local_kmax = grid_data->nzones[2];
  int local_imax_1 = local_imax + 1;
  int local_jmax_1 = local_jmax + 1;

  double *dx = &grid_data->deltas[0][0];
  double *dy = &grid_data->deltas[1][0];
  double * dz = &grid_data->deltas[2][0];

  SubTVec &psi_lf = *gd_set->psi_lf;
  SubTVec &psi_fr = *gd_set->psi_fr;
  SubTVec &psi_bo = *gd_set->psi_bo;

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
  int il = (extent.inc_i > 0 ? 0 : 1);
  int jf = (extent.inc_j > 0 ? 0 : 1);
  int kb = (extent.inc_k > 0 ? 0 : 1);
  int ir = !il;
  int jb = !jf;
  int kt = !kb;

  std::vector<Grid_Sweep_Block> const &idxset =
      grid_data->octant_indexset[octant];
  std::vector<double> xcos_dxi_all(local_imax);
  std::vector<double> ycos_dyj_all(local_jmax);
  std::vector<double> zcos_dzk_all(local_kmax);

  for (int d = 0; d < num_directions; ++d) {
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

    for (int block_idx = 0; block_idx < idxset.size(); ++block_idx) {
      Grid_Sweep_Block const &block = idxset[block_idx];

#ifdef KRIPKE_USE_OPENMP
#pragma omp parallel for
#endif
      for (int group = 0; group < num_groups; ++group) {
        double * KRESTRICT  psi_d_g = gd_set->psi->ptr(group, d, 0);
        double *  KRESTRICT rhs_d_g = gd_set->rhs->ptr(group, d, 0);
        double *  KRESTRICT psi_lf_d_g = psi_lf.ptr(group, d, 0);
        double *  KRESTRICT psi_fr_d_g = psi_fr.ptr(group, d, 0);
        double * KRESTRICT  psi_bo_d_g = psi_bo.ptr(group, d, 0);
        double *  KRESTRICT i_plane_d_g = &i_plane(group, d, 0);
        double *  KRESTRICT j_plane_d_g = &j_plane(group, d, 0);
        double *  KRESTRICT k_plane_d_g = &k_plane(group, d, 0);
        double *  KRESTRICT sigt_g = gd_set->sigt->ptr(group, 0, 0);

        /* Copy the angular fluxes incident upon this subdomain */
        for (int k = block.start_k; k != block.end_k; k += block.inc_k) {
          for (int j = block.start_j; j != block.end_j; j += block.inc_j) {
            /* psi_lf has length (local_imax+1)*local_jmax*local_kmax */
            psi_lf_d_g[Left_INDEX(extent.start_i + il, j, k)] =
                i_plane_d_g[I_PLANE_INDEX(j, k)];
          }
        }

        for (int k = block.start_k; k != block.end_k; k += block.inc_k) {
          for (int i = block.start_i; i != block.end_i; i += block.inc_i) {
            /* psi_fr has length local_imax*(local_jmax+1)*local_kmax */
            psi_fr_d_g[Front_INDEX(i, extent.start_j + jf,
                k)] = j_plane_d_g[J_PLANE_INDEX(i, k)];
          }
        }

        for (int j = block.start_j; j != block.end_j; j += block.inc_j) {
          for (int i = block.start_i; i != block.end_i; i += block.inc_i) {
            /* psi_bo has length local_imax*local_jmax*(local_kmax+1) */
            psi_bo_d_g[Bottom_INDEX(i, j, extent.start_k + kb)] =
                k_plane_d_g[K_PLANE_INDEX(i, j)];
          }
        }


        for (int k = block.start_k; k != block.end_k; k += block.inc_k) {
          double zcos_dzk = zcos_dzk_all[k];
          for (int j = block.start_j; j != block.end_j; j += block.inc_j) {
            double ycos_dyj = ycos_dyj_all[j];
            for (int i = block.start_i; i != block.end_i; i += block.inc_i) {
              double xcos_dxi = xcos_dxi_all[i];

              /* Calculate new zonal flux */
              double psi_d_g_z = (rhs_d_g[Zonal_INDEX(i, j, k)]
                  + psi_lf_d_g[Left_INDEX(i+il, j, k )] * xcos_dxi
                  + psi_fr_d_g[Front_INDEX(i, j+jf, k )] * ycos_dyj
                  + psi_bo_d_g[Bottom_INDEX(i, j, k+kb )] * zcos_dzk)
                  / (xcos_dxi + ycos_dyj + zcos_dzk
                      + sigt_g[Zonal_INDEX(i, j, k)]);
              psi_d_g[Zonal_INDEX(i, j, k)] = psi_d_g_z;
              /* Apply diamond-difference relationships */
              psi_lf_d_g[Left_INDEX(i+ir, j, k )] = 2.0 * psi_d_g_z
                  - psi_lf_d_g[Left_INDEX(i+il, j, k)];
              psi_fr_d_g[Front_INDEX(i, j+jb, k )] = 2.0 * psi_d_g_z
                  - psi_fr_d_g[Front_INDEX(i, j+jf, k )];
              psi_bo_d_g[Bottom_INDEX(i, j, k+kt )] = 2.0 * psi_d_g_z
                  - psi_bo_d_g[Bottom_INDEX(i, j, k+kb)];
            }
          }
        }

        /* Copy the angular fluxes exiting this subdomain */
        for (int k = block.start_k; k != block.end_k; k += block.inc_k) {
          for (int j = block.start_j; j != block.end_j; j += block.inc_j) {
            i_plane_d_g[I_PLANE_INDEX(j, k)] =
                psi_lf_d_g[Left_INDEX(extent.end_i+ir-extent.inc_i, j, k)];
          }
        }
        for (int k = block.start_k; k != block.end_k; k += block.inc_k) {
          for (int i = block.start_i; i != block.end_i; i += block.inc_i) {
            j_plane_d_g[J_PLANE_INDEX(i, k)] =
                psi_fr_d_g[Front_INDEX(i, extent.end_j+jb-extent.inc_j, k)];
          }
        }
        for (int j = block.start_j; j != block.end_j; j += block.inc_j) {
          for (int i = block.start_i; i != block.end_i; i += block.inc_i) {
            k_plane_d_g[K_PLANE_INDEX(i, j)] =
                psi_bo_d_g[Bottom_INDEX(i, j, extent.end_k+kt-extent.inc_k)];
          }
        }

      }
    } // group
  } // direction

}

