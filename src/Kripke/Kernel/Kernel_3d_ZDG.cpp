#include<Kripke/Kernel/Kernel_3d_ZDG.h>
#include<Kripke/User_Data.h>
#include<Kripke/SubTVec.h>
#include<Kripke/LMat.h>

Kernel_3d_ZDG::Kernel_3d_ZDG() {

}

Kernel_3d_ZDG::~Kernel_3d_ZDG() {

}

Nesting_Order Kernel_3d_ZDG::nestingPsi(void) const {
  return NEST_ZDG;
}

Nesting_Order Kernel_3d_ZDG::nestingPhi(void) const {
  return NEST_ZDG;
}

void Kernel_3d_ZDG::scattering(Grid_Data *grid_data) {
  int num_moments = grid_data->num_moments;
  int num_groups = grid_data->phi->groups;
  int num_zones = grid_data->num_zones;

  double ***phi_in = grid_data->phi->data;
  double ***phi_out = grid_data->phi_out->data;

  for (int zone = 0; zone < num_zones; zone++) {

    double **phi_in_z = phi_in[zone];
    double **phi_out_z = phi_out[zone];

    int m0 = 0;
    // Begin loop over scattering moments
    for (int n = 0; n < num_moments; n++) {

      int num_m = grid_data->ell->numM(n);

      for (int m = 0; m < num_m; m++) {

          double * __restrict__ phi_out_z_nm = phi_out_z[m + m0];
          double * __restrict__ phi_in_z_nm = phi_in_z[m + m0];

          // Loop over source group
          for (int g = 0; g < num_groups; g++) {

            // Loop over destination group
            for (int gp = 0; gp < num_groups; gp++) {

            // Evaluate sigs  for this (n,g,gp) triplet
            //evalSigmaS(grid_data, n, g, gp);

            // Get variables
            //double *sig_s = &grid_data->sig_s[0];
              double sig_s = 0;

            phi_out_z_nm[g] += sig_s * phi_in_z_nm[g];
          } //gp

        } // g
      } // m

      m0 += num_m;
    } // n
  } // z
}

void Kernel_3d_ZDG::LTimes(Grid_Data *grid_data) {
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
      for (int z = 0; z < num_zones; z++) {
        double **psi_z = psi[z];
        double **phi_z = phi[z];

        for (int n = 0; n < num_moments; n++) {
          double ** __restrict__ phi_z_n = phi_z + n * n + n;
          double **ell_n = ell[n];
          for (int m = -n; m <= n; m++) {
            double * __restrict__ ell_n_m = ell_n[m + n];
            double * __restrict__ phi_z_nm = phi_z_n[m];

            for (int group = 0; group < num_local_groups; ++group) {
              phi_z_nm[group + group0] = 0;
            }

            for (int d = 0; d < num_local_directions; d++) {
              double * __restrict__ psi_z_d = psi_z[d];
              double ell_n_m_d = ell_n_m[d];

              for (int group = 0; group < num_local_groups; ++group) {
                double psi_z_d_g = psi_z_d[group];
                phi_z_nm[group + group0] += ell_n_m_d * psi_z_d_g;
              }
            }
          }
        }
      }

    } // Direction Set
  } // Group Set
}

void Kernel_3d_ZDG::LPlusTimes(Grid_Data *grid_data) {
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
      for (int z = 0; z < num_zones; z++) {
        double **psi_z = rhs[z];
        double **phi_z = phi_out[z];

        for (int d = 0; d < num_local_directions; d++) {
          double **ell_plus_d = ell_plus[d];
          double * __restrict__ psi_z_d = psi_z[d];

          for (int group = 0; group < num_local_groups; ++group) {
            psi_z_d[group] = 0.0;
          }

          for (int n = 0; n < num_moments; n++) {
            double *ell_plus_d_n = ell_plus_d[n];

            double **phi_z_n = phi_z + n * n;

            for (int m = -n; m <= n; m++) {
              double * __restrict__ phi_z_n_m = phi_z_n[m + n] + group0;
              double ell_plus_d_n_m = ell_plus_d_n[m + n];
              for (int group = 0; group < num_local_groups; ++group) {
                psi_z_d[group] += ell_plus_d_n_m * phi_z_n_m[group];
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

void Kernel_3d_ZDG::sweep(Grid_Data *grid_data, Group_Dir_Set *gd_set,
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

  SubTVec psi_lf(nestingPsi(), num_groups, num_directions,
      (local_imax + 1) * local_jmax * local_kmax);
  SubTVec psi_fr(nestingPsi(), num_groups, num_directions,
      local_imax * (local_jmax + 1) * local_kmax);
  SubTVec psi_bo(nestingPsi(), num_groups, num_directions,
      local_imax * local_jmax * (local_kmax + 1));

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

  /* Copy the angular fluxes incident upon this subdomain */
  for (int k = 0; k < local_kmax; k++) {
    for (int j = 0; j < local_jmax; j++) {
      double ** psi_lf_z = psi_lf.data[Left_INDEX(istartz+il, j, k)];
      double ** i_plane_z = i_plane[I_PLANE_INDEX(j, k)];
      for (int d = 0; d < num_directions; ++d) {
        double * __restrict__ psi_lf_z_d = psi_lf_z[d];
        double * __restrict__ i_plane_z_d = i_plane_z[d];
        for (int group = 0; group < num_groups; ++group) {
          psi_lf_z_d[group] = i_plane_z_d[group];
        }
      }
    }
  }

  for (int k = 0; k < local_kmax; k++) {
    for (int i = 0; i < local_imax; i++) {
      double ** psi_fr_z = psi_fr.data[Front_INDEX(i, jstartz+jf, k)];
      double ** j_plane_z = j_plane[J_PLANE_INDEX(i, k)];

      for (int d = 0; d < num_directions; ++d) {
        double * __restrict__ psi_fr_z_d = psi_fr_z[d];
        double * __restrict__ j_plane_d_z = j_plane_z[d];
        for (int group = 0; group < num_groups; ++group) {
          psi_fr_z_d[group] = j_plane_d_z[group];
        }
      }
    }
  }

  for (int j = 0; j < local_jmax; j++) {
    for (int i = 0; i < local_imax; i++) {
      double ** psi_bo_z = psi_bo.data[Bottom_INDEX(i, j, kstartz+ kb)];
      double ** k_plane_z = k_plane[K_PLANE_INDEX(i, j)];

      for (int d = 0; d < num_directions; ++d) {
        double * __restrict__ psi_bo_z_d = psi_bo_z[d];
        double * __restrict__ k_plane_z_d = k_plane_z[d];
        for (int group = 0; group < num_groups; ++group) {
          psi_bo_z_d[group] = k_plane_z_d[group];
        }
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
        double **psi_z = psi[z];
        double **rhs_z = rhs[z];

        double **psi_lf_zil = psi_lf.data[Left_INDEX(i+il, j, k)];
        double **psi_lf_zir = psi_lf.data[Left_INDEX(i+ir, j, k)];

        double **psi_fr_zjf = psi_fr.data[Front_INDEX(i, j+jf, k)];
        double **psi_fr_zjb = psi_fr.data[Front_INDEX(i, j+jb, k)];

        double **psi_bo_zkb = psi_bo.data[Bottom_INDEX(i, j, k+kb)];
        double **psi_bo_zkt = psi_bo.data[Bottom_INDEX(i, j, k+kt)];

        double **psi_internal_all_z = psi_internal_all[z];
        double **i_plane_z = i_plane[I_PLANE_INDEX(j, k)];
        double **j_plane_z = j_plane[J_PLANE_INDEX(i, k)];
        double **k_plane_z = k_plane[K_PLANE_INDEX(i, j)];

        for (int d = 0; d < num_directions; ++d) {
          double xcos = direction[d].xcos;
          double ycos = direction[d].ycos;
          double zcos = direction[d].zcos;

          double zcos_dzk = 2.0 * zcos / dzk;
          double ycos_dyj = 2.0 * ycos / dyj;
          double xcos_dxi = 2.0 * xcos / dxi;

          double * __restrict__ psi_z_d = psi_z[d];
          double * __restrict__ rhs_z_d = rhs_z[d];

          double * __restrict__ psi_lf_zil_d = psi_lf_zil[d];
          double * __restrict__ psi_lf_zir_d = psi_lf_zir[d];

          double * __restrict__ psi_fr_zjf_d = psi_fr_zjf[d];
          double * __restrict__ psi_fr_zjb_d = psi_fr_zjb[d];

          double * __restrict__ psi_bo_zkb_d = psi_bo_zkb[d];
          double * __restrict__ psi_bo_zkt_d = psi_bo_zkt[d];

          double * __restrict__ psi_internal_all_z_d = psi_internal_all_z[d];
          double * __restrict__ i_plane_z_d = i_plane_z[d];
          double * __restrict__ j_plane_z_d = j_plane_z[d];
          double * __restrict__ k_plane_z_d = k_plane_z[d];

          double * __restrict__ sigt_z = sigt[z];

          double *psi_int_lf = psi_internal_all_z_d;
          double *psi_int_fr = psi_internal_all_z_d;
          double *psi_int_bo = psi_internal_all_z_d;

          for (int group = 0; group < num_groups; ++group) {

            /* Add internal surface source data */
            psi_lf_zil_d[group] += psi_int_lf[group];
            psi_fr_zjf_d[group] += psi_int_fr[group];
            psi_bo_zkb_d[group] += psi_int_bo[group];

            /* Calculate new zonal flux */
            double psi_z_d_g = (rhs_z_d[group] + psi_lf_zil_d[group] * xcos_dxi
                + psi_fr_zjf_d[group] * ycos_dyj
                + psi_bo_zkb_d[group] * zcos_dzk)
                / (xcos_dxi + ycos_dyj + zcos_dzk + sigt_z[group]);

            psi_z_d[group] = psi_z_d_g;

            /* Apply diamond-difference relationships */
            psi_lf_zir_d[group] = 2.0 * psi_z_d_g - psi_lf_zil_d[group];
            psi_fr_zjb_d[group] = 2.0 * psi_z_d_g - psi_fr_zjf_d[group];
            psi_bo_zkt_d[group] = 2.0 * psi_z_d_g - psi_bo_zkb_d[group];
          }
        }
      }
    }

  }

  /* Copy the angular fluxes exiting this subdomain */
  for (int k = 0; k < local_kmax; k++) {
    for (int j = 0; j < local_jmax; j++) {
      double ** psi_lf_z = psi_lf.data[Left_INDEX(istopz+ir, j, k)];
      double ** i_plane_z = i_plane[I_PLANE_INDEX(j, k)];
      for (int d = 0; d < num_directions; ++d) {
        double * __restrict__ psi_lf_z_d = psi_lf_z[d];
        double * __restrict__ i_plane_z_d = i_plane_z[d];
        for (int group = 0; group < num_groups; ++group) {
          i_plane_z_d[group] = psi_lf_z_d[group];
        }
      }
    }
  }

  for (int k = 0; k < local_kmax; k++) {
    for (int i = 0; i < local_imax; i++) {
      double ** psi_fr_z = psi_fr.data[Front_INDEX(i, jstopz+jb, k)];
      double ** j_plane_z = j_plane[J_PLANE_INDEX(i, k)];

      for (int d = 0; d < num_directions; ++d) {
        double * __restrict__ psi_fr_z_d = psi_fr_z[d];
        double * __restrict__ j_plane_d_z = j_plane_z[d];
        for (int group = 0; group < num_groups; ++group) {
          j_plane_d_z[group] = psi_fr_z_d[group];
        }
      }
    }
  }

  for (int j = 0; j < local_jmax; j++) {
    for (int i = 0; i < local_imax; i++) {
      double ** psi_bo_z = psi_bo.data[Bottom_INDEX(i, j, kstopz+kt)];
      double ** k_plane_z = k_plane[K_PLANE_INDEX(i, j)];

      for (int d = 0; d < num_directions; ++d) {
        double * __restrict__ psi_bo_z_d = psi_bo_z[d];
        double * __restrict__ k_plane_z_d = k_plane_z[d];
        for (int group = 0; group < num_groups; ++group) {
          k_plane_z_d[group] = psi_bo_z_d[group];
        }
      }
    }
  }
}

