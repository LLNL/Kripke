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

void Kernel_3d_DZG::scattering(Grid_Data *grid_data) {
  int num_moments = grid_data->num_moments;
  int num_groups = grid_data->phi->groups;
  int num_zones = grid_data->num_zones;

  double ***phi_in = grid_data->phi->data;
  double ***phi_out = grid_data->phi_out->data;

  double sig_s = 0;

  // Begin loop over scattering moments
  int m0 = 0;
  for (int n = 0; n < num_moments; n++) {
    int num_m = grid_data->ell->numM(n);

    for (int m = 0; m < num_m; m++) {

      double **phi_in_nm = phi_in[m0 + m];
      double **phi_out_nm = phi_out[m0 + m];

      for (int zone = 0; zone < num_zones; zone++) {

        double * __restrict__ phi_out_nm_z = phi_out_nm[zone];
        double * __restrict__ phi_in_nm_z = phi_in_nm[zone];

        // Loop over destination group
        for (int gp = 0; gp < num_groups; gp++) {

          // Loop over source group

          for (int g = 0; g < num_groups; g++) {

            // Evaluate sigs  for this (n,g,gp) triplet
            //evalSigmaS(grid_data, n, g, gp);

            // Get variables
            //double *sig_s = &grid_data->sig_s[0];

            phi_out_nm_z[g] += sig_s * phi_in_nm_z[g];
          } // g

        } // gp
      } // z

    } // m

    m0 += num_m;
  } // n
}

void Kernel_3d_DZG::LTimes(Grid_Data *grid_data) {
  // Outer parameters
  double ***phi = grid_data->phi->data;
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

      // Get Variables
      double ***psi = gd_set.psi->data;

      /* 3D Cartesian Geometry */
      for (int n = 0; n < num_moments; n++) {
        double ***phi_n = phi + n * n;
        double **ell_n = ell[n];

#ifdef KRIPKE_USE_OPENMP
#pragma omp parallel for
#endif
        for (int m = -n; m <= n; m++) {
          double **phi_n_m = phi_n[m + n];
          double *ell_n_m = ell_n[m + n];
          for (int d = 0; d < num_local_directions; d++) {
            double **psi_d = psi[d];
            double ell_n_m_d = ell_n_m[d + dir0];
            for (int z = 0; z < num_zones; z++) {
              double * __restrict__ psi_d_z = psi_d[z];
              double * __restrict__ phi_n_z = phi_n_m[z];
              for (int group = 0; group < num_local_groups; ++group) {
                double psi_d_g_z = psi_d_z[group];
                phi_n_z[group + group0] += ell_n_m_d * psi_d_g_z;
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
  double ***phi_out = grid_data->phi_out->data;
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
      double ***rhs = gd_set.rhs->data;
      gd_set.rhs->clear(0.0);

      /* 3D Cartesian Geometry */
      for (int d = 0; d < num_local_directions; d++) {
        double **rhs_d = rhs[d];
        double **ell_plus_d = ell_plus[d + dir0];

        for (int n = 0; n < num_moments; n++) {
          double ***phi_out_n = phi_out + n * n;
          double *ell_plus_d_n = ell_plus_d[n];

          for (int m = 0; m <= 2 * n; m++) {
            double **phi_out_n_m = phi_out_n[m];
            double ell_plus_d_n_m = ell_plus_d_n[m];

#ifdef KRIPKE_USE_OPENMP
#pragma omp parallel for
#endif
            for (int z = 0; z < num_zones; z++) {
              double * __restrict__ rhs_d_z = rhs_d[z];
              double * __restrict__ phi_n_m_z = phi_out_n_m[z] + group0;

              for (int group = 0; group < num_local_groups; ++group) {
                rhs_d_z[group] += ell_plus_d_n_m * phi_n_m_z[group];
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

  double * __restrict__ dx = &grid_data->deltas[0][0];
  double * __restrict__ dy = &grid_data->deltas[1][0];
  double * __restrict__ dz = &grid_data->deltas[2][0];

  SubTVec &psi_lf = *gd_set->psi_lf;
  SubTVec &psi_fr = *gd_set->psi_fr;
  SubTVec &psi_bo = *gd_set->psi_bo;

  // Alias the MPI data with a SubTVec for the face data
  SubTVec i_plane_v(nestingPsi(), num_groups, num_directions,
      local_jmax * local_kmax, i_plane_ptr);
  SubTVec j_plane_v(nestingPsi(), num_groups, num_directions,
      local_imax * local_kmax, j_plane_ptr);
  SubTVec k_plane_v(nestingPsi(), num_groups, num_directions,
      local_imax * local_jmax, k_plane_ptr);
  double ***i_plane = i_plane_v.data;
  double ***j_plane = j_plane_v.data;
  double ***k_plane = k_plane_v.data;

  double ***psi_internal_all = gd_set->psi_internal->data;

  double ***psi = gd_set->psi->data;
  double ***rhs = gd_set->rhs->data;
  double **sigt = gd_set->sigt->data[0];

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
  for (int d = 0; d < num_directions; ++d) {
    double **psi_d = psi[d];
    double **rhs_d = rhs[d];
    double **psi_lf_d = psi_lf.data[d];
    double **psi_fr_d = psi_fr.data[d];
    double **psi_bo_d = psi_bo.data[d];
    double **psi_internal_all_d = psi_internal_all[d];
    double **i_plane_d = i_plane[d];
    double **j_plane_d = j_plane[d];
    double **k_plane_d = k_plane[d];

    double xcos = direction[d].xcos;
    double ycos = direction[d].ycos;
    double zcos = direction[d].zcos;

    /* Copy the angular fluxes incident upon this subdomain */
    for (int k = 0; k < local_kmax; k++) {
      for (int j = 0; j < local_jmax; j++) {
        /* psi_lf has length (local_imax+1)*local_jmax*local_kmax */
        double * __restrict__ psi_lf_d_z =
            psi_lf_d[Left_INDEX(extent.start_i+il, j, k)];
        double * __restrict__ i_plane_d_z = i_plane_d[I_PLANE_INDEX(j, k)];
        for (int group = 0; group < num_groups; ++group) {
          psi_lf_d_z[group] = i_plane_d_z[group];
        }
      }
    }

    for (int k = 0; k < local_kmax; k++) {
      for (int i = 0; i < local_imax; i++) {
        /* psi_fr has length local_imax*(local_jmax+1)*local_kmax */
        double * __restrict__ psi_fr_d_z =
            psi_fr_d[Front_INDEX(i, extent.start_j+jf, k)];
        double * __restrict__ j_plane_d_z = j_plane_d[J_PLANE_INDEX(i, k)];
        for (int group = 0; group < num_groups; ++group) {
          psi_fr_d_z[group] = j_plane_d_z[group];
        }
      }
    }

    for (int j = 0; j < local_jmax; j++) {
      for (int i = 0; i < local_imax; i++) {
        double * __restrict__ psi_bo_d_z =
            psi_bo_d[Bottom_INDEX(i, j, extent.start_k+ kb)];
        double * __restrict__ k_plane_d_z = k_plane_d[K_PLANE_INDEX(i, j)];
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
            double * __restrict__ psi_d_z = psi_d[z];
            double * __restrict__ rhs_d_z = rhs_d[z];

            double * __restrict__ psi_lf_d_zil =
                psi_lf_d[Left_INDEX(i+il, j, k)];
            double * __restrict__ psi_lf_d_zir =
                psi_lf_d[Left_INDEX(i+ir, j, k)];

            double * __restrict__ psi_fr_d_zjf =
                psi_fr_d[Front_INDEX(i, j+jf, k)];
            double * __restrict__ psi_fr_d_zjb =
                psi_fr_d[Front_INDEX(i, j+jb, k)];

            double * __restrict__ psi_bo_d_zkb =
                psi_bo_d[Bottom_INDEX(i, j, k+kb)];
            double * __restrict__ psi_bo_d_zkt =
                psi_bo_d[Bottom_INDEX(i, j, k+kt)];

            double * __restrict__ psi_internal_all_d_z = psi_internal_all_d[z];
            double * __restrict__ i_plane_d_z = i_plane_d[I_PLANE_INDEX(j, k)];
            double * __restrict__ j_plane_d_z = j_plane_d[J_PLANE_INDEX(i, k)];
            double * __restrict__ k_plane_d_z = k_plane_d[K_PLANE_INDEX(i, j)];

            double * __restrict__ sigt_z = sigt[z];

            for (int group = 0; group < num_groups; ++group) {
              double *psi_int_lf = psi_internal_all_d_z;
              double *psi_int_fr = psi_internal_all_d_z;
              double *psi_int_bo = psi_internal_all_d_z;

              /* Add internal surface source data */
              psi_lf_d_zil[group] += psi_int_lf[group];
              psi_fr_d_zjf[group] += psi_int_fr[group];
              psi_bo_d_zkb[group] += psi_int_bo[group];

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
        double * __restrict__ psi_lf_d_z = psi_lf_d[Left_INDEX(extent.end_i+ir-extent.inc_i, j, k)];
        double * __restrict__ i_plane_d_z = i_plane_d[I_PLANE_INDEX(j, k)];
        for (int group = 0; group < num_groups; ++group) {
          i_plane_d_z[group] = psi_lf_d_z[group];
        }
      }
    }

    for (int k = 0; k < local_kmax; k++) {
      for (int i = 0; i < local_imax; i++) {
        double * __restrict__ psi_fr_d_z =
            psi_fr_d[Front_INDEX(i, extent.end_j+jb-extent.inc_j, k)];
        double * __restrict__ j_plane_d_z = j_plane_d[J_PLANE_INDEX(i, k)];
        for (int group = 0; group < num_groups; ++group) {
          j_plane_d_z[group] = psi_fr_d_z[group];
        }
      }
    }

    for (int j = 0; j < local_jmax; j++) {
      for (int i = 0; i < local_imax; i++) {
        double * __restrict__ psi_bo_d_z =
            psi_bo_d[Bottom_INDEX(i, j, extent.end_k+kt-extent.inc_k)];
        double * __restrict__ k_plane_d_z = k_plane_d[K_PLANE_INDEX(i, j)];
        for (int group = 0; group < num_groups; ++group) {
          k_plane_d_z[group] = psi_bo_d_z[group];
        }
      }
    }

  }

}

