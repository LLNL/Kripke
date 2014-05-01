#include "../Param.h"


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

void zdg_nmd_SweepDD_3d(zdg_nmd_Param &p) {
  int num_directions = p.num_directions;
  int num_zones = p.num_zones;
  int num_groups = p.num_groups;
  Direction *direction = p.direction;

  int local_imax = p.nzones[0];
  int local_jmax = p.nzones[1];
  int local_kmax = p.nzones[2];
  int local_imax_1 = local_imax + 1;
  int local_jmax_1 = local_jmax + 1;

  double * __restrict__ dx = p.deltas[0];
  double * __restrict__ dy = p.deltas[1];
  double * __restrict__ dz = p.deltas[2];

  zdg_TVec psi_lf_v(num_groups, num_directions,
                    (local_imax+1)*local_jmax*local_kmax);
  zdg_TVec psi_fr_v(num_groups, num_directions,
                    local_imax*(local_jmax+1)*local_kmax);
  zdg_TVec psi_bo_v(num_groups, num_directions,
                    local_imax*local_jmax*(local_kmax+1));
  double ***psi_lf = psi_lf_v.data;
  double ***psi_fr = psi_fr_v.data;
  double ***psi_bo = psi_bo_v.data;

  zdg_TVec i_plane_v(num_groups, num_directions, local_jmax*local_kmax);
  zdg_TVec j_plane_v(num_groups, num_directions, local_imax*local_kmax);
  zdg_TVec k_plane_v(num_groups, num_directions, local_imax*local_jmax);
  double ***i_plane = i_plane_v.data;
  double ***j_plane = j_plane_v.data;
  double ***k_plane = k_plane_v.data;

  zdg_TVec psi_internal(num_groups, num_directions, num_zones);
  double ***psi_internal_all = psi_internal.data;

  double ***psi = p.psi.data;
  double ***rhs = p.rhs.data;
  double **sigt = p.sigt.data;

  // All directions have same id,jd,kd, since we are modeling an "Angle Set"
  // So pull that information out now
  int istartz, istopz, in, il, ir;
  int id = direction[0].id;
  int jd = direction[0].jd;
  int kd = direction[0].kd;
  if(id > 0){
    istartz = 0; istopz = local_imax-1; in = 1; il = 0; ir = 1;
  }
  else {
    istartz = local_imax-1; istopz = 0; in = -1; il = 1; ir = 0;
  }

  int jstartz, jstopz, jn, jf, jb;
  if(jd > 0){
    jstartz = 0; jstopz = local_jmax-1; jn = 1; jf = 0; jb = 1;
  }
  else {
    jstartz = local_jmax-1; jstopz = 0; jn = -1; jf = 1; jb = 0;
  }

  int kstartz, kstopz, kn, kb, kt;
  if(kd > 0){
    kstartz = 0; kstopz = local_kmax-1; kn =  1; kb = 0; kt = 1;
  }
  else {
    kstartz = local_kmax-1; kstopz = 0; kn = -1; kb = 1; kt = 0;
  }


  /* Copy the angular fluxes incident upon this subdomain */
  for(int k=0; k<local_kmax; k++){
    for(int j=0; j<local_jmax; j++){
      double ** psi_lf_z =
        psi_lf[Left_INDEX(istartz+il, j, k)];
      double ** i_plane_z =
        i_plane[I_PLANE_INDEX(j, k)];
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for(int d = 0; d < num_directions; ++d){
        double * __restrict__ psi_lf_z_d = psi_lf_z[d];
        double * __restrict__ i_plane_z_d = i_plane_z[d];
        for(int group = 0; group < num_groups; ++group){
          psi_lf_z_d[group] = i_plane_z_d[group];
        }
      }
    }
  }

  for(int k=0; k<local_kmax; k++){
    for(int i=0; i<local_imax; i++){
      double ** psi_fr_z =
        psi_fr[Front_INDEX(i, jstartz+jf, k)];
      double ** j_plane_z =
        j_plane[J_PLANE_INDEX(i, k)];

#ifdef _OPENMP
#pragma omp parallel for
#endif
      for(int d = 0; d < num_directions; ++d){
        double * __restrict__ psi_fr_z_d = psi_fr_z[d];
        double * __restrict__ j_plane_d_z = j_plane_z[d];
        for(int group = 0; group < num_groups; ++group){
          psi_fr_z_d[group] = j_plane_d_z[group];
        }
      }
    }
  }

  for(int j=0; j<local_jmax; j++){
    for(int i=0; i<local_imax; i++){
      double ** psi_bo_z =
        psi_bo[Bottom_INDEX(i, j, kstartz+ kb)];
      double ** k_plane_z =
        k_plane[K_PLANE_INDEX(i, j)];

#ifdef _OPENMP
#pragma omp parallel for
#endif
      for(int d = 0; d < num_directions; ++d){
        double * __restrict__ psi_bo_z_d =psi_bo_z[d];
        double * __restrict__ k_plane_z_d = k_plane_z[d];
        for(int group = 0; group < num_groups; ++group){
          psi_bo_z_d[group] = k_plane_z_d[group];
        }
      }
    }
  }

  /*  Perform transport sweep of the grid 1 cell at a time.   */
  for(int k=kstartz; std::abs(k-kstartz)<local_kmax; k+=kn){
    double dzk = dz[k+1];
    for(int j=jstartz; std::abs(j-jstartz)<local_jmax; j+=jn){
      double dyj = dy[j+1];
      for(int i=istartz; std::abs(i-istartz)<local_imax; i+=in){
        double dxi = dx[i+1];

        int z = Zonal_INDEX(i, j, k);
        double **psi_z = psi[z];
        double **rhs_z = rhs[z];

        double **psi_lf_zil = psi_lf[Left_INDEX(i+il, j, k)];
        double **psi_lf_zir = psi_lf[Left_INDEX(i+ir, j, k)];

        double **psi_fr_zjf = psi_fr[Front_INDEX(i, j+jf, k)];
        double **psi_fr_zjb = psi_fr[Front_INDEX(i, j+jb, k)];

        double **psi_bo_zkb = psi_bo[Bottom_INDEX(i, j, k+kb)];
        double **psi_bo_zkt = psi_bo[Bottom_INDEX(i, j, k+kt)];

        double **psi_internal_all_z = psi_internal_all[z];
        double **i_plane_z = i_plane[I_PLANE_INDEX(j, k)];
        double **j_plane_z = j_plane[J_PLANE_INDEX(i, k)];
        double **k_plane_z = k_plane[K_PLANE_INDEX(i, j)];

#ifdef _OPENMP
#pragma omp parallel for
#endif
        for(int d = 0; d < num_directions; ++d){
          double xcos = direction[d].xcos;
          double ycos = direction[d].ycos;
          double zcos = direction[d].zcos;

          double zcos_dzk = 2.0*zcos/dzk;
          double ycos_dyj = 2.0*ycos/dyj;
          double xcos_dxi = 2.0*xcos/dxi;

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

          for(int group = 0; group < num_groups; ++group){

            /* Add internal surface source data */
            psi_lf_zil_d[group] += psi_int_lf[group];
            psi_fr_zjf_d[group] += psi_int_fr[group];
            psi_bo_zkb_d[group] += psi_int_bo[group];

            /* Calculate new zonal flux */
            double psi_z_d_g =
              (rhs_z_d[group]
               + psi_lf_zil_d[group]*xcos_dxi
               + psi_fr_zjf_d[group]*ycos_dyj
               + psi_bo_zkb_d[group]*zcos_dzk)/
              (xcos_dxi + ycos_dyj + zcos_dzk + sigt_z[group]);

            psi_z_d[group] = psi_z_d_g;

            /* Apply diamond-difference relationships */
            psi_lf_zir_d[group] = 2.0*psi_z_d_g - psi_lf_zil_d[group];
            psi_fr_zjb_d[group] = 2.0*psi_z_d_g - psi_fr_zjf_d[group];
            psi_bo_zkt_d[group] = 2.0*psi_z_d_g - psi_bo_zkb_d[group];
          }
        }
      }
    }

  }

  /* Copy the angular fluxes exiting this subdomain */
  for(int k=0; k<local_kmax; k++){
    for(int j=0; j<local_jmax; j++){
      double ** psi_lf_z =
        psi_lf[Left_INDEX(istopz+ir, j, k)];
      double ** i_plane_z =
        i_plane[I_PLANE_INDEX(j, k)];
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for(int d = 0; d < num_directions; ++d){
        double * __restrict__ psi_lf_z_d = psi_lf_z[d];
        double * __restrict__ i_plane_z_d = i_plane_z[d];
        for(int group = 0; group < num_groups; ++group){
          i_plane_z_d[group] = psi_lf_z_d[group];
        }
      }
    }
  }

  for(int k=0; k<local_kmax; k++){
    for(int i=0; i<local_imax; i++){
      double ** psi_fr_z =
        psi_fr[Front_INDEX(i, jstopz+jb, k)];
      double ** j_plane_z =
        j_plane[J_PLANE_INDEX(i, k)];

#ifdef _OPENMP
#pragma omp parallel for
#endif
      for(int d = 0; d < num_directions; ++d){
        double * __restrict__ psi_fr_z_d = psi_fr_z[d];
        double * __restrict__ j_plane_d_z = j_plane_z[d];
        for(int group = 0; group < num_groups; ++group){
          j_plane_d_z[group] = psi_fr_z_d[group];
        }
      }
    }
  }

  for(int j=0; j<local_jmax; j++){
    for(int i=0; i<local_imax; i++){
      double ** psi_bo_z =
        psi_bo[Bottom_INDEX(i, j, kstopz+kt)];
      double ** k_plane_z =
        k_plane[K_PLANE_INDEX(i, j)];

#ifdef _OPENMP
#pragma omp parallel for
#endif
      for(int d = 0; d < num_directions; ++d){
        double * __restrict__ psi_bo_z_d =psi_bo_z[d];
        double * __restrict__ k_plane_z_d = k_plane_z[d];
        for(int group = 0; group < num_groups; ++group){
          k_plane_z_d[group] = psi_bo_z_d[group];
        }
      }
    }
  }


}
