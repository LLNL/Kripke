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

  double ***phi_in = grid_data->phi->data;
  double ***phi_out = grid_data->phi_out->data;

  double sig_s = 0;

  // Loop over destination group
  for (int gp = 0; gp < num_groups; gp++) {

    for (int g = 0; g < num_groups; g++) {

      double **phi_in_g = phi_in[g];
      double **phi_out_g = phi_out[g];

      for (int zone = 0; zone < num_zones; zone++) {
        int m0 = 0;

        double * __restrict__ phi_out_g_z = phi_out_g[zone];
        double * __restrict__ phi_in_g_z = phi_in_g[zone];

        // Begin loop over scattering moments
        for (int n = 0; n < num_moments; n++) {

          int num_m = grid_data->ell->numM(n);

          // Loop over source group

          // Evaluate sigs  for this (n,g,gp) triplet
          //evalSigmaS(grid_data, n, g, gp);

          // Get variables
          //double *sig_s = &grid_data->sig_s[0];

          for (int m = 0; m < num_m; m++) {
            phi_out_g_z[m + m0] += sig_s * phi_in_g_z[m + m0];
          } //m

          m0 += num_m;

        } // n
      } // z
    } // g
  } // gp
}

void Kernel_3d_GZD::LTimes(Grid_Data *grid_data) {
  // Outer parameters
  double ***phi = grid_data->phi->data;
  double ***ell = grid_data->ell->data;
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
      double ***psi = gd_set.psi->data;

      /* 3D Cartesian Geometry */
      for (int group = 0; group < num_local_groups; ++group) {
        double **psi_zonal = psi[group];
        double **phi_g = phi[group + group0];

        for (int z = 0; z < num_zones; z++) {
          double * __restrict__ psi_g_z = psi_zonal[z];
          double * __restrict__ phi_g_z = phi_g[z];

          for (int n = 0; n < num_moments; n++) {
            double * __restrict__ phi_g_z_n = phi_g_z + (n * n);
            double **ell_n = ell[n];
            int nn = 2 * n;
            for (int m = 0; m <= nn; m++) {

              double* __restrict__ ell_n_m = ell_n[m];
              double phi_g_z_nm = 0.0;

              for (int d = 0; d < num_local_directions; d++) {
                double ell_n_m_d = ell_n_m[d + dir0];
                double psi_g_z_d = psi_g_z[d];
                phi_g_z_nm += ell_n_m_d * psi_g_z_d;
              }

              phi_g_z_n[m] = phi_g_z_nm;
            }
          }
        }
      }

    } // Direction Set
  } // Group Set
}

void Kernel_3d_GZD::LPlusTimes(Grid_Data *grid_data) {
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

      /* 3D Cartesian Geometry */
      for (int group = 0; group < num_local_groups; ++group) {
        double **phi_out_g = phi_out[group + group0];
        double **rhs_g = rhs[group];

        for (int i = 0; i < num_zones; i++) {
          double const * __restrict__ phi_out_g_z = phi_out_g[i];
          double * __restrict__ rhs_g_z = rhs_g[i];
          for (int d = 0; d < num_local_directions; d++) {
            double **ell_plus_d = ell_plus[d + dir0];
            double psi_g_z_d = 0.0;

            for (int n = 0; n < num_moments; n++) {
              int nn = n * n;
              int n2 = 2 * n;
              double const * __restrict__ ell_plus_d_n = ell_plus_d[n];
              double const * __restrict__ phi_g_z_n = phi_out_g_z + nn;

              double psi_g_z_d_m = 0.0;
              for (int m = 0; m <= n2; m++) {
                psi_g_z_d_m += ell_plus_d_n[m] * phi_g_z_n[m];
              }
              psi_g_z_d += psi_g_z_d_m;
            }
            rhs_g_z[d] = psi_g_z_d;
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

  SubTVec psi_internal(nestingPsi(), num_groups, num_directions, num_zones);
  double ***psi_internal_all = psi_internal.data;

  double ***psi = gd_set->psi->data;
  double ***rhs = gd_set->rhs->data;
  double **sigt = gd_set->sigt->data[0];

  // All directions have same id,jd,kd, since these are all one Direction Set
  // So pull that information out now
  int istartz, istopz, in, il, ir;
  int id = direction[0].id;
  int jd = direction[0].jd;
  int kd = direction[0].kd;
  if (id > 0) {
    istartz = 0;
    istopz = local_imax - 1;
    in = 1;
    il = 0;
    ir = 1;
  } else {
    istartz = local_imax - 1;
    istopz = 0;
    in = -1;
    il = 1;
    ir = 0;
  }

  int jstartz, jstopz, jn, jf, jb;
  if (jd > 0) {
    jstartz = 0;
    jstopz = local_jmax - 1;
    jn = 1;
    jf = 0;
    jb = 1;
  } else {
    jstartz = local_jmax - 1;
    jstopz = 0;
    jn = -1;
    jf = 1;
    jb = 0;
  }

  int kstartz, kstopz, kn, kb, kt;
  if (kd > 0) {
    kstartz = 0;
    kstopz = local_kmax - 1;
    kn = 1;
    kb = 0;
    kt = 1;
  } else {
    kstartz = local_kmax - 1;
    kstopz = 0;
    kn = -1;
    kb = 1;
    kt = 0;
  }

  for (int group = 0; group < num_groups; ++group) {
    double **psi_g = psi[group];
    double **rhs_g = rhs[group];
    double **psi_lf_g = psi_lf.data[group];
    double **psi_fr_g = psi_fr.data[group];
    double **psi_bo_g = psi_bo.data[group];
    double **psi_internal_all_g = psi_internal_all[group];
    double **i_plane_g = i_plane[group];
    double **j_plane_g = j_plane[group];
    double **k_plane_g = k_plane[group];
    double * __restrict__ sigt_g = sigt[group];

    /* Copy the angular fluxes incident upon this subdomain */
    for (int k = 0; k < local_kmax; k++) {
      for (int j = 0; j < local_jmax; j++) {
        /* psi_lf has length (local_imax+1)*local_jmax*local_kmax */
        double * __restrict__ psi_lf_g_z =
            psi_lf_g[Left_INDEX(istartz+il, j, k)];
        double * __restrict__ i_plane_g_z = i_plane_g[I_PLANE_INDEX(j, k)];
        for (int d = 0; d < num_directions; ++d) {
          psi_lf_g_z[d] = i_plane_g_z[d];
        }
      }
    }

    for (int k = 0; k < local_kmax; k++) {
      for (int i = 0; i < local_imax; i++) {
        /* psi_fr has length local_imax*(local_jmax+1)*local_kmax */
        double * __restrict__ psi_fr_g_z =
            psi_fr_g[Front_INDEX(i, jstartz+jf, k)];
        double * __restrict__ j_plane_g_z = j_plane_g[J_PLANE_INDEX(i, k)];
        for (int d = 0; d < num_directions; ++d) {
          psi_fr_g_z[d] = j_plane_g_z[d];
        }
      }
    }

    for (int j = 0; j < local_jmax; j++) {
      for (int i = 0; i < local_imax; i++) {
        double * __restrict__ psi_bo_g_z =
            psi_bo_g[Bottom_INDEX(i, j, kstartz+ kb)];
        double * __restrict__ k_plane_g_z = k_plane_g[K_PLANE_INDEX(i, j)];
        for (int d = 0; d < num_directions; ++d) {
          psi_bo_g_z[d] = k_plane_g_z[d];
        }
      }
    }

    /*  Perform transport sweep of the grid 1 cell at a time.   */
    for (int k = kstartz; std::abs(k - kstartz) < local_kmax; k += kn) {
      double dzk = dz[k + 1];
      for (int j = jstartz; std::abs(j - jstartz) < local_jmax; j += jn) {
        double dyj = dy[j + 1];
        for (int i = istartz; std::abs(i - istartz) < local_imax; i += in) {
          double dxi = dx[i + 1];

          int z = Zonal_INDEX(i, j, k);
          double * __restrict__ psi_g_z = psi_g[z];
          double * __restrict__ rhs_g_z = rhs_g[z];

          double * __restrict__ psi_lf_g_zil = psi_lf_g[Left_INDEX(i+il, j, k)];
          double * __restrict__ psi_lf_g_zir = psi_lf_g[Left_INDEX(i+ir, j, k)];

          double * __restrict__ psi_fr_g_zjf = psi_fr_g[Front_INDEX(i, j+jf, k)];
          double * __restrict__ psi_fr_g_zjb = psi_fr_g[Front_INDEX(i, j+jb, k)];

          double * __restrict__ psi_bo_g_zkb =
              psi_bo_g[Bottom_INDEX(i, j, k+kb)];
          double * __restrict__ psi_bo_g_zkt =
              psi_bo_g[Bottom_INDEX(i, j, k+kt)];

          double * __restrict__ psi_internal_all_g_z = psi_internal_all_g[z];
          double * __restrict__ i_plane_g_z = i_plane_g[I_PLANE_INDEX(j, k)];
          double * __restrict__ j_plane_g_z = j_plane_g[J_PLANE_INDEX(i, k)];
          double * __restrict__ k_plane_g_z = k_plane_g[K_PLANE_INDEX(i, j)];

          for (int d = 0; d < num_directions; ++d) {
            double xcos = direction[d].xcos;
            double ycos = direction[d].ycos;
            double zcos = direction[d].zcos;

            double zcos_dzk = 2.0 * zcos / dzk;
            double ycos_dyj = 2.0 * ycos / dyj;
            double xcos_dxi = 2.0 * xcos / dxi;

            double *psi_int_lf = psi_internal_all_g_z;
            double *psi_int_fr = psi_internal_all_g_z;
            double *psi_int_bo = psi_internal_all_g_z;

            /* Add internal surface source data */
            psi_lf_g_zil[d] += psi_int_lf[d];
            psi_fr_g_zjf[d] += psi_int_fr[d];
            psi_bo_g_zkb[d] += psi_int_bo[d];

            /* Calculate new zonal flux */
            double psi_g_z_d =
                (rhs_g_z[d] + psi_lf_g_zil[d] * xcos_dxi
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

    /* Copy the angular fluxes exiting this subdomain */
    for (int k = 0; k < local_kmax; k++) {
      for (int j = 0; j < local_jmax; j++) {
        double * __restrict__ psi_lf_g_z = psi_lf_g[Left_INDEX(istopz+ir, j, k)];
        double * __restrict__ i_plane_g_z = i_plane_g[I_PLANE_INDEX(j, k)];
        for (int d = 0; d < num_directions; ++d) {
          psi_lf_g_z[d] = i_plane_g_z[d];
        }
      }
    }

    for (int k = 0; k < local_kmax; k++) {
      for (int i = 0; i < local_imax; i++) {
        double * __restrict__ psi_fr_g_z =
            psi_fr_g[Front_INDEX(i, jstopz+jb, k)];
        double * __restrict__ j_plane_g_z = j_plane_g[J_PLANE_INDEX(i, k)];
        for (int d = 0; d < num_directions; ++d) {
          psi_fr_g_z[d] = j_plane_g_z[d];
        }
      }
    }

    for (int j = 0; j < local_jmax; j++) {
      for (int i = 0; i < local_imax; i++) {
        double * __restrict__ psi_bo_g_z =
            psi_bo_g[Bottom_INDEX(i, j, kstopz+kt)];
        double * __restrict__ k_plane_g_z = k_plane_g[K_PLANE_INDEX(i, j)];
        for (int d = 0; d < num_directions; ++d) {
          k_plane_g_z[d] = psi_bo_g_z[d];
        }
      }
    }

  } // group

}

