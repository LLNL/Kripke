#include<Kripke/Kernel/Kernel_3d_DZG.h>
#include<Kripke/User_Data.h>
#include<Kripke/SubTVec.h>
#include<Kripke/LMat.h>

Kernel_3d_DZG::Kernel_3d_DZG() {

}

Kernel_3d_DZG::~Kernel_3d_DZG() {

}

Nesting_Order Kernel_3d_DZG::nestingPsi(void) const {
  return NEST_DZG;
}

Nesting_Order Kernel_3d_DZG::nestingPhi(void) const {
  return NEST_DZG;
}

void Kernel_3d_DZG::LTimes(Grid_Data *grid_data) {
  // Outer parameters
  double ***ell = grid_data->ell->data;
  int num_zones = grid_data->num_zones;
  int num_moments = grid_data->num_moments;

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
      for (int n = 0; n < num_moments; n++) {
        double **ell_n = ell[n];
        for (int m = -n; m <= n; m++) {
          int nm_offset = n*n + n + m;
          double *ell_n_m = ell_n[m + n];
          for (int d = 0; d < num_local_directions; d++) {
            double ell_n_m_d = ell_n_m[d + dir0];
            for (int z = 0; z < num_zones; z++) {

              double *phi = grid_data->phi->ptr(group0, nm_offset, z);
              double *psi = gd_set.psi->ptr(0, d, z);
              for (int group = 0; group < num_local_groups; ++group) {
                phi[group] += ell_n_m_d * psi[group];
              }

            }
          }
        }
      }

    } // Direction Set
  } // Group Set
}

void Kernel_3d_DZG::LPlusTimes(Grid_Data *grid_data) {
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
      for (int d = 0; d < num_local_directions; d++) {
        double **ell_plus_d = ell_plus[d + dir0];

        for (int n = 0; n < num_moments; n++) {
          double *ell_plus_d_n = ell_plus_d[n];

          for (int m = -n; m <= n; m++) {
            int nm_offset = n*n + n + m;
            double ell_plus_d_n_m = ell_plus_d_n[n+m];

            for (int z = 0; z < num_zones; z++) {
              double *phi_out = grid_data->phi_out->ptr(group0, nm_offset, z);
              double *rhs = gd_set.rhs->ptr(0, d, z);

              for (int group = 0; group < num_local_groups; ++group) {
                rhs[group] += ell_plus_d_n_m * phi_out[group];
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

void Kernel_3d_DZG::sweep(Grid_Data *grid_data, Group_Dir_Set *gd_set,
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

  for (int d = 0; d < num_directions; ++d) {
    double xcos = direction[d].xcos;
    double ycos = direction[d].ycos;
    double zcos = direction[d].zcos;

    /* Copy the angular fluxes incident upon this subdomain */
    for (int k = 0; k < local_kmax; k++) {
      for (int j = 0; j < local_jmax; j++) {
        /* psi_lf has length (local_imax+1)*local_jmax*local_kmax */
        double * psi_lf_d_z = psi_lf.ptr(0, d, Left_INDEX(extent.start_i+il, j, k));
        double * i_plane_d_z = i_plane.ptr(0, d, I_PLANE_INDEX(j, k));
        for (int group = 0; group < num_groups; ++group) {
          psi_lf_d_z[group] = i_plane_d_z[group];
        }
      }
    }

    for (int k = 0; k < local_kmax; k++) {
      for (int i = 0; i < local_imax; i++) {
        /* psi_fr has length local_imax*(local_jmax+1)*local_kmax */
        double * psi_fr_d_z = psi_fr.ptr(0, d, Front_INDEX(i, extent.start_j+jf, k));
        double * j_plane_d_z = j_plane.ptr(0, d, J_PLANE_INDEX(i, k));
        for (int group = 0; group < num_groups; ++group) {
          psi_fr_d_z[group] = j_plane_d_z[group];
        }
      }
    }

    for (int j = 0; j < local_jmax; j++) {
      for (int i = 0; i < local_imax; i++) {
        double * psi_bo_d_z = psi_bo.ptr(0, d, Bottom_INDEX(i, j, extent.start_k+ kb));
        double * k_plane_d_z = k_plane.ptr(0, d, K_PLANE_INDEX(i, j));
        for (int group = 0; group < num_groups; ++group) {
          psi_bo_d_z[group] = k_plane_d_z[group];
        }
      }
    }

    /*  Perform transport sweep of the grid 1 cell at a time.   */
    for (int block_idx = 0; block_idx < idxset.size(); ++block_idx) {
      Grid_Sweep_Block const &block = idxset[block_idx];
      for (int k = block.start_k; k != block.end_k; k += block.inc_k) {
        double dzk = dz[k + 1];
        double zcos_dzk = 2.0 * zcos / dzk;
        for (int j = block.start_j; j != block.end_j; j += block.inc_j) {
          double dyj = dy[j + 1];
          double ycos_dyj = 2.0 * ycos / dyj;
          for (int i = block.start_i; i != block.end_i; i += block.inc_i) {
            double dxi = dx[i + 1];
            double xcos_dxi = 2.0 * xcos / dxi;

            int z = Zonal_INDEX(i, j, k);
            double * psi_d_z = gd_set->psi->ptr(0, d, z);
            double * rhs_d_z = gd_set->rhs->ptr(0, d, z);

            double * psi_lf_d_zil = psi_lf.ptr(0, d, Left_INDEX(i+il, j, k));
            double * psi_lf_d_zir = psi_lf.ptr(0, d, Left_INDEX(i+ir, j, k));

            double * psi_fr_d_zjf = psi_fr.ptr(0, d, Front_INDEX(i, j+jf, k));
            double * psi_fr_d_zjb = psi_fr.ptr(0, d, Front_INDEX(i, j+jb, k));

            double * psi_bo_d_zkb = psi_bo.ptr(0, d, Bottom_INDEX(i, j, k+kb));
            double * psi_bo_d_zkt = psi_bo.ptr(0, d, Bottom_INDEX(i, j, k+kt));

            double * i_plane_d_z = i_plane.ptr(0, d, I_PLANE_INDEX(j, k));
            double * j_plane_d_z = j_plane.ptr(0, d, J_PLANE_INDEX(i, k));
            double * k_plane_d_z = k_plane.ptr(0, d, K_PLANE_INDEX(i, j));

            double * sigt_z = gd_set->sigt->ptr(0, 0, z);

            for (int group = 0; group < num_groups; ++group) {
              /* Calculate new zonal flux */
              double psi_d_z_g = (rhs_d_z[group]
                  + psi_lf_d_zil[group] * xcos_dxi
                  + psi_fr_d_zjf[group] * ycos_dyj
                  + psi_bo_d_zkb[group] * zcos_dzk)
                  / (xcos_dxi + ycos_dyj + zcos_dzk + sigt_z[group]);

              psi_d_z[group] = psi_d_z_g;

              /* Apply diamond-difference relationships */
              psi_lf_d_zir[group] = 2.0 * psi_d_z_g - psi_lf_d_zil[group];
              psi_fr_d_zjb[group] = 2.0 * psi_d_z_g - psi_fr_d_zjf[group];
              psi_bo_d_zkt[group] = 2.0 * psi_d_z_g - psi_bo_d_zkb[group];
            }
          }
        }
      }
    }

    /* Copy the angular fluxes exiting this subdomain */
    for (int k = 0; k < local_kmax; k++) {
      for (int j = 0; j < local_jmax; j++) {
        double * psi_lf_d_z =
            psi_lf.ptr(0, d, Left_INDEX(extent.end_i+ir-extent.inc_i, j, k));
        double * i_plane_d_z = i_plane.ptr(0, d, I_PLANE_INDEX(j, k));
        for (int group = 0; group < num_groups; ++group) {
          i_plane_d_z[group] = psi_lf_d_z[group];
        }
      }
    }

    for (int k = 0; k < local_kmax; k++) {
      for (int i = 0; i < local_imax; i++) {
        double * psi_fr_d_z =
            psi_fr.ptr(0, d, Front_INDEX(i, extent.end_j+jb-extent.inc_j, k));
        double * j_plane_d_z = j_plane.ptr(0, d, J_PLANE_INDEX(i, k));
        for (int group = 0; group < num_groups; ++group) {
          j_plane_d_z[group] = psi_fr_d_z[group];
        }
      }
    }

    for (int j = 0; j < local_jmax; j++) {
      for (int i = 0; i < local_imax; i++) {
        double * psi_bo_d_z =
            psi_bo.ptr(0, d, Bottom_INDEX(i, j, extent.end_k+kt-extent.inc_k));
        double * k_plane_d_z = k_plane.ptr(0, d, K_PLANE_INDEX(i, j));
        for (int group = 0; group < num_groups; ++group) {
          k_plane_d_z[group] = psi_bo_d_z[group];
        }
      }
    }

  }

}

