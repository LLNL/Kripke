#include<Kripke/Kernel/Kernel_3d_GZD.h>
#include<Kripke/User_Data.h>
#include<Kripke/SubTVec.h>
#include<Kripke/LMat.h>

Kernel_3d_GZD::Kernel_3d_GZD() {

}

Kernel_3d_GZD::~Kernel_3d_GZD() {

}

Nesting_Order Kernel_3d_GZD::nestingPsi(void) const {
  return NEST_GZD;
}

Nesting_Order Kernel_3d_GZD::nestingPhi(void) const {
  return NEST_GZD;
}

void Kernel_3d_GZD::scattering(Grid_Data *grid_data) {
  int num_moments = grid_data->num_moments;
  int num_groups = grid_data->phi->groups;
  int num_zones = grid_data->num_zones;

  double sig_s = 0;

  // Loop over destination group
  for (int gp = 0; gp < num_groups; gp++) {
    for (int g = 0; g < num_groups; g++) {
      for (int zone = 0; zone < num_zones; zone++) {

        // Begin loop over scattering moments
        double *phi = grid_data->phi->ptr(g, 0, zone);
        double *phi_out = grid_data->phi_out->ptr(g, 0, zone);
        int nm_offset = 0;
        for (int n = 0; n < num_moments; n++) {

          int num_m = grid_data->ell->numM(n);
          // Loop over source group

          // Evaluate sigs  for this (n,g,gp) triplet
          //evalSigmaS(grid_data, n, g, gp);

          // Get variables
          //double *sig_s = &grid_data->sig_s[0];

          for (int m = 0; m < num_m; m++) {
            phi_out[nm_offset + m] += sig_s * phi[nm_offset + m];
          } //m

          nm_offset += num_m;
        } // n
      } // z
    } // g
  } // gp
}

void Kernel_3d_GZD::LTimes(Grid_Data *grid_data) {
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
#ifdef KRIPKE_USE_OPENMP
#pragma omp parallel for
#endif
        for (int z = 0; z < num_zones; z++) {
          double *psi = gd_set.psi->ptr(group, 0, z);

          int nm_offset = 0;
          double *phi = grid_data->phi->ptr(group0+group, 0, z);
          for (int n = 0; n < num_moments; n++) {
            double **ell_n = ell[n];

            for (int m = -n; m <= n; m++) {
              double *ell_n_m = ell_n[n+m] + dir0;
              double phi_acc = 0.0;

              for (int d = 0; d < num_local_directions; d++) {
                double ell_n_m_d = ell_n_m[d];
                double psi_g_z_d = psi[d];
                phi_acc += ell_n_m_d * psi_g_z_d;
              }

              phi[nm_offset] += phi_acc;
              nm_offset ++;
            }
          }
        }
      }

    } // Direction Set
  } // Group Set
}

void Kernel_3d_GZD::LPlusTimes(Grid_Data *grid_data) {
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

      /* 3D Cartesian Geometry */
      for (int group = 0; group < num_local_groups; ++group) {
#ifdef KRIPKE_USE_OPENMP
#pragma omp parallel for
#endif
        for (int z = 0; z < num_zones; z++) {
          double *rhs = gd_set.rhs->ptr(group, 0, z);
          double const *phi_out = grid_data->phi_out->ptr(group0+group, 0, z);
          for (int d = 0; d < num_local_directions; d++) {
            double **ell_plus_d = ell_plus[d + dir0];

            double psi_acc = 0.0;
            int nm_offset = 0;
            for (int n = 0; n < num_moments; n++) {
              double const * ell_plus_d_n = ell_plus_d[n];
              for (int m = -n; m <= n; m++) {
                psi_acc += ell_plus_d_n[n+m] * phi_out[nm_offset+n+m];
              }
              nm_offset += 2*n+1;
            }
            rhs[d] = psi_acc;
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

void Kernel_3d_GZD::sweep(Grid_Data *grid_data, Group_Dir_Set *gd_set,
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

#ifdef KRIPKE_USE_OPENMP
#pragma omp parallel for
#endif
  for (int group = 0; group < num_groups; ++group) {
    double *sigt_g = gd_set->sigt->ptr(group,0 , 0);

    /* Copy the angular fluxes incident upon this subdomain */
    for (int k = 0; k < local_kmax; k++) {
      for (int j = 0; j < local_jmax; j++) {
        /* psi_lf has length (local_imax+1)*local_jmax*local_kmax */
        double * psi_lf_g_z = psi_lf.ptr(group, 0, Left_INDEX(extent.start_i+il, j, k));
        double * i_plane_g_z = i_plane.ptr(group, 0, I_PLANE_INDEX(j, k));
        for (int d = 0; d < num_directions; ++d) {
          psi_lf_g_z[d] = i_plane_g_z[d];
        }
      }
    }

    for (int k = 0; k < local_kmax; k++) {
      for (int i = 0; i < local_imax; i++) {
        /* psi_fr has length local_imax*(local_jmax+1)*local_kmax */
        double * psi_fr_g_z = psi_fr.ptr(group, 0, Front_INDEX(i, extent.start_j+jf, k));
        double * j_plane_g_z = j_plane.ptr(group, 0, J_PLANE_INDEX(i, k));
        for (int d = 0; d < num_directions; ++d) {
          psi_fr_g_z[d] = j_plane_g_z[d];
        }
      }
    }

    for (int j = 0; j < local_jmax; j++) {
      for (int i = 0; i < local_imax; i++) {
        double * psi_bo_g_z = psi_bo.ptr(group, 0, Bottom_INDEX(i, j, extent.start_k+ kb));
        double * k_plane_g_z = k_plane.ptr(group, 0, K_PLANE_INDEX(i, j));
        for (int d = 0; d < num_directions; ++d) {
          psi_bo_g_z[d] = k_plane_g_z[d];
        }
      }
    }

    /*  Perform transport sweep of the grid 1 cell at a time.   */
    for (int block_idx = 0; block_idx < idxset.size(); ++block_idx) {
      Grid_Sweep_Block const &block = idxset[block_idx];
      for (int k = block.start_k; k != block.end_k; k += block.inc_k) {
        double dzk = dz[k + 1];
        double two_dz = 2.0 / dzk;
        for (int j = block.start_j; j != block.end_j; j += block.inc_j) {
          double dyj = dy[j + 1];
          double two_dy = 2.0 / dyj;
          for (int i = block.start_i; i != block.end_i; i += block.inc_i) {
            double dxi = dx[i + 1];
            double two_dx = 2.0 / dxi;

            int z = Zonal_INDEX(i, j, k);
            double *psi_g_z = gd_set->psi->ptr(group, 0, z);
            double *rhs_g_z = gd_set->rhs->ptr(group, 0, z);

            double * psi_lf_g_zil = psi_lf.ptr(group, 0, Left_INDEX(i+il, j, k));
            double * psi_lf_g_zir = psi_lf.ptr(group, 0, Left_INDEX(i+ir, j, k));

            double * psi_fr_g_zjf = psi_fr.ptr(group, 0, Front_INDEX(i, j+jf, k));
            double * psi_fr_g_zjb = psi_fr.ptr(group, 0, Front_INDEX(i, j+jb, k));

            double * psi_bo_g_zkb = psi_bo.ptr(group, 0, Bottom_INDEX(i, j, k+kb));
            double * psi_bo_g_zkt = psi_bo.ptr(group, 0, Bottom_INDEX(i, j, k+kt));

            double *psi_internal_all_g_z = gd_set->psi_internal->ptr(group, 0, z);
            double * i_plane_g_z = i_plane.ptr(group, 0, I_PLANE_INDEX(j, k));
            double * j_plane_g_z = j_plane.ptr(group, 0, J_PLANE_INDEX(i, k));
            double * k_plane_g_z = k_plane.ptr(group, 0, K_PLANE_INDEX(i, j));

            for (int d = 0; d < num_directions; ++d) {
              double xcos = direction[d].xcos;
              double ycos = direction[d].ycos;
              double zcos = direction[d].zcos;

              double zcos_dzk = zcos * two_dz;
              double ycos_dyj = ycos * two_dy;
              double xcos_dxi = xcos * two_dx;

              double *psi_int_lf = psi_internal_all_g_z;
              double *psi_int_fr = psi_internal_all_g_z;
              double *psi_int_bo = psi_internal_all_g_z;

              /* Add internal surface source data */
              psi_lf_g_zil[d] += psi_int_lf[d];
              psi_fr_g_zjf[d] += psi_int_fr[d];
              psi_bo_g_zkb[d] += psi_int_bo[d];

              /* Calculate new zonal flux */
              double psi_g_z_d = (rhs_g_z[d] + psi_lf_g_zil[d] * xcos_dxi
                  + psi_fr_g_zjf[d] * ycos_dyj + psi_bo_g_zkb[d] * zcos_dzk)
                  / (xcos_dxi + ycos_dyj + zcos_dzk
                      + sigt_g[Zonal_INDEX(i, j, k)]);

              psi_g_z[d] = psi_g_z_d;

              /* Apply diamond-difference relationships */
              psi_lf_g_zir[d] = 2.0 * psi_g_z_d - psi_lf_g_zil[d];
              psi_fr_g_zjb[d] = 2.0 * psi_g_z_d - psi_fr_g_zjf[d];
              psi_bo_g_zkt[d] = 2.0 * psi_g_z_d - psi_bo_g_zkb[d];
            }
          }
        }
      }
    }
    /* Copy the angular fluxes exiting this subdomain */
    for (int k = 0; k < local_kmax; k++) {
      for (int j = 0; j < local_jmax; j++) {
        double * psi_lf_g_z =
            psi_lf.ptr(group, 0, Left_INDEX(extent.end_i-extent.inc_i+ir, j, k));
        double * i_plane_g_z = i_plane.ptr(group, 0, I_PLANE_INDEX(j, k));
        for (int d = 0; d < num_directions; ++d) {
          i_plane_g_z[d] = psi_lf_g_z[d];
        }
      }
    }

    for (int k = 0; k < local_kmax; k++) {
      for (int i = 0; i < local_imax; i++) {
        double * psi_fr_g_z =
            psi_fr.ptr(group, 0, Front_INDEX(i, extent.end_j-extent.inc_j+jb, k));
        double * j_plane_g_z = j_plane.ptr(group, 0, J_PLANE_INDEX(i, k));
        for (int d = 0; d < num_directions; ++d) {
          j_plane_g_z[d] = psi_fr_g_z[d];
        }
      }
    }

    for (int j = 0; j < local_jmax; j++) {
      for (int i = 0; i < local_imax; i++) {
        double * psi_bo_g_z =
            psi_bo.ptr(group, 0, Bottom_INDEX(i, j, extent.end_k-extent.inc_k+kt));
        double * k_plane_g_z = k_plane.ptr(group, 0, K_PLANE_INDEX(i, j));
        for (int d = 0; d < num_directions; ++d) {
          k_plane_g_z[d] = psi_bo_g_z[d];
        }
      }
    }

  } // group

}

