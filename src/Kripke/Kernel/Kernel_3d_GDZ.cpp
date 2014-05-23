#include<Kripke/Kernel/Kernel_3d_GDZ.h>
#include<Kripke/User_Data.h>
#include<Kripke/SubTVec.h>
#include<Kripke/LMat.h>

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

void Kernel_3d_GDZ::scattering(Grid_Data *grid_data) {
  int num_moments = grid_data->num_moments;
  int num_groups = grid_data->phi->groups;
  int num_zones = grid_data->num_zones;

    // Loop over destination group
  for (int gp = 0; gp < num_groups; gp++) {

    // Begin loop over scattering moments
    int nm_offset0 = 0;
    for (int n = 0; n < num_moments; n++) {

      int num_m = grid_data->ell->numM(n);

      // Loop over source group
      for (int g = 0; g < num_groups; g++) {

        // Evaluate sigs  for this (n,g,gp) triplet
        evalSigmaS(grid_data, n, g, gp);

        // Get variables
        double *sig_s = &grid_data->sig_s[0];

        int nm_offset = nm_offset0;
        for (int m = 0; m < num_m; m++) {
          double *phi = grid_data->phi->ptr(g, nm_offset, 0);
          double *phi_out = grid_data->phi_out->ptr(g, nm_offset, 0);

          for (int zone = 0; zone < num_zones; zone++) {
            phi_out[zone] += sig_s[zone] * phi[zone];
          }

          nm_offset ++;
        } // m
      } // g
      nm_offset0 += num_m;
    } // n
  } // gp
}

void Kernel_3d_GDZ::LTimes(Grid_Data *grid_data) {
  // Outer parameters
  double ***ell = grid_data->ell->data;
  int num_zones = grid_data->num_zones;
  int num_moments = grid_data->num_moments;

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
      for (int group = 0; group < num_local_groups; ++group) {
        for (int n = 0; n < num_moments; n++) {
          double **ell_n = ell[n];
          for (int m = -n; m <= n; m++) {
            int nm_offset = n*n + n + m;
            double * ell_n_m = ell_n[m + n];

            double * KRESTRICT phi = grid_data->phi->ptr(group0+group, nm_offset, 0);
            double * KRESTRICT psi = gd_set.psi->ptr(group, 0, 0);

            for (int d = 0; d < num_local_directions; d++) {
              double ell_n_m_d = ell_n_m[d + dir0];
              for (int i = 0; i < num_zones; i++) {
                phi[i] += ell_n_m_d * psi[i];
              }
              psi += num_zones;
            }

          }
        }
      }

    } // Direction Set
  } // Group Set
}

void Kernel_3d_GDZ::LPlusTimes(Grid_Data *grid_data) {
  // Outer parameters
  double ***ell_plus = grid_data->ell_plus->data;
  int num_zones = grid_data->num_zones;
  int num_moments = grid_data->num_moments;

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
      for (int group = 0; group < num_local_groups; ++group) {
        for (int d = 0; d < num_local_directions; d++) {
          double **ell_plus_d = ell_plus[d + dir0];

          for (int n = 0; n < num_moments; n++) {
            double *ell_plus_d_n = ell_plus_d[n];
            for (int m = -n; m <= n; m++) {
              int nm_offset = n*n + n + m;

              double ell_plus_d_n_m = ell_plus_d_n[m + n];

              double *phi_out = grid_data->phi_out->ptr(group0+group, nm_offset, 0);
              double *rhs = gd_set.rhs->ptr(group, d, 0);

              for (int z = 0; z < num_zones; z++) {
                rhs[z] += ell_plus_d_n_m * phi_out[z];
              }
            }
          }
        }
      }

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

void Kernel_3d_GDZ::sweep(Grid_Data *grid_data, Group_Dir_Set *gd_set,
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

  double * dx = &grid_data->deltas[0][0];
  double * dy = &grid_data->deltas[1][0];
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

  for (int group = 0; group < num_groups; ++group) {
    double * sigt_g = gd_set->sigt->ptr(group, 0, 0);

    for (int d = 0; d < num_directions; ++d) {
      double * KRESTRICT psi_g_d = gd_set->psi->ptr(group, d, 0);
      double * KRESTRICT rhs_g_d = gd_set->rhs->ptr(group, d, 0);
      double * KRESTRICT psi_lf_g_d = gd_set->psi_lf->ptr(group, d, 0);
      double * KRESTRICT psi_fr_g_d = gd_set->psi_fr->ptr(group, d, 0);
      double * KRESTRICT psi_bo_g_d = gd_set->psi_bo->ptr(group, d, 0);
      double * KRESTRICT psi_internal_all_g_d = gd_set->psi_internal->ptr(group, d, 0);
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

      double * KRESTRICT psi_int_lf = psi_internal_all_g_d;
      double * KRESTRICT psi_int_fr = psi_internal_all_g_d;
      double * KRESTRICT psi_int_bo = psi_internal_all_g_d;

      /* Copy the angular fluxes incident upon this subdomain */
      for (int k = 0; k < local_kmax; k++) {
        for (int j = 0; j < local_jmax; j++) {
          /* psi_lf has length (local_imax+1)*local_jmax*local_kmax */
          psi_lf_g_d[Left_INDEX(extent.start_i+il, j,
              k)] = i_plane_g_d[I_PLANE_INDEX(j, k)];
        }
      }

      for (int k = 0; k < local_kmax; k++) {
        for (int i = 0; i < local_imax; i++) {
          /* psi_fr has length local_imax*(local_jmax+1)*local_kmax */
          psi_fr_g_d[Front_INDEX(i, extent.start_j+jf,
              k)] = j_plane_g_d[J_PLANE_INDEX(i, k)];
        }
      }

      for (int j = 0; j < local_jmax; j++) {
        for (int i = 0; i < local_imax; i++) {
          /* psi_bo has length local_imax*local_jmax*(local_kmax+1) */
          psi_bo_g_d[Bottom_INDEX(i, j, extent.start_k+
              kb)] = k_plane_g_d[K_PLANE_INDEX(i, j)];
        }
      }

      /*  Perform transport sweep of the grid 1 cell at a time.   */
      for (int block_idx = 0; block_idx < idxset.size(); ++block_idx) {
        Grid_Sweep_Block const &block = idxset[block_idx];
        for (int k = block.start_k; k != block.end_k; k += block.inc_k) {
          double zcos_dzk = zcos_dzk_all[k];
          for (int j = block.start_j; j != block.end_j; j += block.inc_j) {
            double ycos_dyj = ycos_dyj_all[j];
            for (int i = block.start_i; i != block.end_i; i += block.inc_i) {
              double xcos_dxi = xcos_dxi_all[i];

              /* Add internal surface source data */
              psi_lf_g_d[Left_INDEX(i+il, j,
                  k)] += psi_int_lf[Zonal_INDEX(i, j, k)];
              psi_fr_g_d[Front_INDEX(i, j+jf,
                  k)] += psi_int_fr[Zonal_INDEX(i, j, k)];
              psi_bo_g_d[Bottom_INDEX(i, j, k+
                  kb)] += psi_int_bo[Zonal_INDEX(i, j, k)];

              /* Calculate new zonal flux */
              double psi_g_d_z = (rhs_g_d[Zonal_INDEX(i, j, k)]
                  + psi_lf_g_d[Left_INDEX(i+il, j, k )] * xcos_dxi
                  + psi_fr_g_d[Front_INDEX(i, j+jf, k )] * ycos_dyj
                  + psi_bo_g_d[Bottom_INDEX(i, j, k+kb )] * zcos_dzk)
                  / (xcos_dxi + ycos_dyj + zcos_dzk
                      + sigt_g[Zonal_INDEX(i, j, k)]);
              psi_g_d[Zonal_INDEX(i, j, k)] = psi_g_d_z;
              /* Apply diamond-difference relationships */
              psi_lf_g_d[Left_INDEX(i+ir, j, k )] = 2.0 * psi_g_d_z
                  - psi_lf_g_d[Left_INDEX(i+il,j, k )];
              psi_fr_g_d[Front_INDEX(i, j+jb, k )] = 2.0 * psi_g_d_z
                  - psi_fr_g_d[Front_INDEX(i,j+ jf,k )];
              psi_bo_g_d[Bottom_INDEX(i, j, k+kt )] = 2.0 * psi_g_d_z
                  - psi_bo_g_d[Bottom_INDEX(i, j, k + kb )];
            }
          }
        }

        /* Copy the angular fluxes exiting this subdomain */
        for (int k = 0; k < local_kmax; k++) {
          for (int j = 0; j < local_jmax; j++) {
            i_plane_g_d[I_PLANE_INDEX(j,k)] =
                psi_lf_g_d[Left_INDEX(extent.end_i-extent.inc_i+ir, j, k)];
          }
        }

        for (int k = 0; k < local_kmax; k++) {
          for (int i = 0; i < local_imax; i++) {
            j_plane_g_d[J_PLANE_INDEX(i, k)] =
                psi_fr_g_d[Front_INDEX(i, extent.end_j-extent.inc_j+jb, k)];
          }
        }

        for (int j = 0; j < local_jmax; j++) {
          for (int i = 0; i < local_imax; i++) {
            k_plane_g_d[K_PLANE_INDEX(i, j)] =
                psi_bo_g_d[Bottom_INDEX(i, j, extent.end_k-extent.inc_k+kt)];
          }
        }
      }
    }
  }

}

