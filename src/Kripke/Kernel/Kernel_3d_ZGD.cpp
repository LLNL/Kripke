#include<Kripke/Kernel/Kernel_3d_ZGD.h>
#include<Kripke/User_Data.h>
#include<Kripke/SubTVec.h>
#include<Kripke/LMat.h>

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
      for (int z = 0; z < num_zones; z++) {
        for (int group = 0; group < num_local_groups; ++group) {
          double *phi = grid_data->phi->ptr(group+group0, 0, z);
          double *psi = gd_set.psi->ptr(group, 0, z);

          for (int n = 0; n < num_moments; n++) {
            double **ell_n = ell[n];

            for (int m = -n; m <= n; m++) {
              int nm_offset = n*n + n + m;
              double * ell_n_m = ell[n][m + n];

              double phi_acc = 0.0;
              for (int d = 0; d < num_local_directions; d++) {
                double ell_n_m_d = ell_n_m[d + dir0];
                double psi_z_g_d = psi[d];
                phi_acc += ell_n_m_d * psi_z_g_d;
              }

              *phi += phi_acc;
              phi++;
            }
          }
        }
      }

    } // Direction Set
  } // Group Set
}

void Kernel_3d_ZGD::LPlusTimes(Grid_Data *grid_data) {
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
      for (int z = 0; z < num_zones; z++) {
        for (int group = 0; group < num_local_groups; ++group) {
          double *rhs = gd_set.rhs->ptr(group, 0, z);

          for (int d = 0; d < num_local_directions; d++) {
            double **ell_plus_d = ell_plus[d + dir0];
            double psi_z_g_d = 0.0;

            double *phi_out = grid_data->phi_out->ptr(group+group0, 0, z);

            for (int n = 0; n < num_moments; n++) {
              double * ell_plus_d_n = ell_plus_d[n];

              int mmax = 2 * n + 1;
              for (int m = 0; m < mmax; m++) {
                psi_z_g_d += ell_plus_d_n[m] * phi_out[m];
              }
              phi_out += mmax;
            }

            rhs[d] = psi_z_g_d;
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

void Kernel_3d_ZGD::sweep(Grid_Data *grid_data, Group_Dir_Set *gd_set,
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

  /* Copy the angular fluxes incident upon this subdomain */
  for (int k = 0; k < local_kmax; k++) {
    for (int j = 0; j < local_jmax; j++) {
      double *psi_lf_z_d = psi_lf.ptr(0, 0, Left_INDEX(extent.start_i+il, j, k));
      double *i_plane_z_d = i_plane.ptr(0, 0, I_PLANE_INDEX(j, k));
      for(int off = 0;off < num_directions*num_groups;++ off){
        psi_lf_z_d[off] = i_plane_z_d[off];
      }
    }
  }

  for (int k = 0; k < local_kmax; k++) {
    for (int i = 0; i < local_imax; i++) {
      double *psi_fr_z_d = psi_fr.ptr(0, 0, Front_INDEX(i, extent.start_j+jf, k));
      double *j_plane_z_d = j_plane.ptr(0, 0, J_PLANE_INDEX(i, k));
      for(int off = 0;off < num_directions*num_groups;++ off){
        psi_fr_z_d[off] = j_plane_z_d[off];
      }
    }
  }

  for (int j = 0; j < local_jmax; j++) {
    for (int i = 0; i < local_imax; i++) {
      double *psi_bo_z_d = psi_bo.ptr(0, 0, Bottom_INDEX(i, j, extent.start_k+ kb));
      double *k_plane_z_d = k_plane.ptr(0, 0, K_PLANE_INDEX(i, j));
      for(int off = 0;off < num_directions*num_groups;++ off){
        psi_bo_z_d[off] = k_plane_z_d[off];
      }
    }
  }

  /*  Perform transport sweep of the grid 1 cell at a time.   */
  for (int block_idx = 0; block_idx < idxset.size(); ++block_idx) {
    Grid_Sweep_Block const &block = idxset[block_idx];

    for (int k = block.start_k; k != block.end_k; k += block.inc_k) {
      double dzk = dz[k + 1];
      for (int j = block.start_j; j != block.end_j; j += block.inc_j) {
        double dyj = dy[j + 1];
        for (int i = block.start_i; i != block.end_i; i += block.inc_i) {
          double dxi = dx[i + 1];

          int z = Zonal_INDEX(i, j, k);
          double * sigt_z = gd_set->sigt->ptr(0, 0, z);

          for (int group = 0; group < num_groups; ++group) {

            double * psi_z_g = gd_set->psi->ptr(group, 0, z);
            double * rhs_z_g = gd_set->rhs->ptr(group, 0, z);

            double *psi_lf_zil_g = psi_lf.ptr(group, 0, Left_INDEX(i+il, j, k));
            double *psi_lf_zir_g = psi_lf.ptr(group, 0, Left_INDEX(i+ir, j, k));

            double *psi_fr_zjf_g = psi_fr.ptr(group, 0, Front_INDEX(i, j+jf, k));
            double *psi_fr_zjb_g = psi_fr.ptr(group, 0, Front_INDEX(i, j+jb, k));

            double *psi_bo_zkb_g = psi_bo.ptr(group, 0, Bottom_INDEX(i, j, k+kb));
            double *psi_bo_zkt_g = psi_bo.ptr(group, 0, Bottom_INDEX(i, j, k+kt));

            double * psi_internal_all_z_g = gd_set->psi_internal->ptr(group, 0, z);
            double * i_plane_z_g = i_plane.ptr(group, 0, z);
            double * j_plane_z_g = j_plane.ptr(group, 0, z);
            double * k_plane_z_g = k_plane.ptr(group, 0, z);

            double *psi_int_lf = psi_internal_all_z_g;
            double *psi_int_fr = psi_internal_all_z_g;
            double *psi_int_bo = psi_internal_all_z_g;

            for (int d = 0; d < num_directions; ++d) {

              double xcos = direction[d].xcos;
              double ycos = direction[d].ycos;
              double zcos = direction[d].zcos;

              double zcos_dzk = 2.0 * zcos / dzk;
              double ycos_dyj = 2.0 * ycos / dyj;
              double xcos_dxi = 2.0 * xcos / dxi;

              /* Add internal surface source data */
              psi_lf_zil_g[d] += psi_int_lf[d];
              psi_fr_zjf_g[d] += psi_int_fr[d];
              psi_bo_zkb_g[d] += psi_int_bo[d];

              /* Calculate new zonal flux */
              double psi_z_g_d = (rhs_z_g[d] + psi_lf_zil_g[d] * xcos_dxi
                  + psi_fr_zjf_g[d] * ycos_dyj + psi_bo_zkb_g[d] * zcos_dzk)
                  / (xcos_dxi + ycos_dyj + zcos_dzk + sigt_z[group]);

              psi_z_g[d] = psi_z_g_d;

              /* Apply diamond-difference relationships */
              psi_lf_zir_g[d] = 2.0 * psi_z_g_d - psi_lf_zil_g[d];
              psi_fr_zjb_g[d] = 2.0 * psi_z_g_d - psi_fr_zjf_g[d];
              psi_bo_zkt_g[d] = 2.0 * psi_z_g_d - psi_bo_zkb_g[d];
            }
          }
        }
      }
    }
  }

  /* Copy the angular fluxes exiting this subdomain */
  for (int k = 0; k < local_kmax; k++) {
    for (int j = 0; j < local_jmax; j++) {
      double *psi_lf_z_d = psi_lf.ptr(0, 0, Left_INDEX(extent.start_i+il, j, k));
      double *i_plane_z_d = i_plane.ptr(0, 0, I_PLANE_INDEX(j, k));
      for(int off = 0;off < num_directions*num_groups;++ off){
        i_plane_z_d[off] = psi_lf_z_d[off];
      }
    }
  }

  for (int k = 0; k < local_kmax; k++) {
    for (int i = 0; i < local_imax; i++) {
      double *psi_fr_z_d = psi_fr.ptr(0, 0, Front_INDEX(i, extent.start_j+jf, k));
      double *j_plane_z_d = j_plane.ptr(0, 0, J_PLANE_INDEX(i, k));
      for(int off = 0;off < num_directions*num_groups;++ off){
        j_plane_z_d[off] = psi_fr_z_d[off];
      }
    }
  }

  for (int j = 0; j < local_jmax; j++) {
    for (int i = 0; i < local_imax; i++) {
      double *psi_bo_z_d = psi_bo.ptr(0, 0, Bottom_INDEX(i, j, extent.start_k+ kb));
      double *k_plane_z_d = k_plane.ptr(0, 0, K_PLANE_INDEX(i, j));
      for(int off = 0;off < num_directions*num_groups;++ off){
        k_plane_z_d[off] = psi_bo_z_d[off];
      }
    }
  }
}

