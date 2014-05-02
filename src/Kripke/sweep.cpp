#include "transport_headers.h"

#include<cstdlib>

/* Sweep routine for Diamond-Difference */
/* Macros for offsets with fluxes on cell faces */
#define Zonal_INDEX(i, j, k) (i) + (local_imax)*(j) \
  + (local_imax)*(local_jmax)*(k)
#define Left_INDEX(i, j, k) (i) + (local_imax_1)*(j) \
  + (local_imax_1)*(local_jmax)*(k)
#define Front_INDEX(i, j, k) (i) + (local_imax)*(j) \
  + (local_imax)*(local_jmax_1)*(k)
#define Bottom_INDEX(i, j, k) (i) + (local_imax)*(j) \
  + (local_imax)*(local_jmax)*(k)
#define I_PLANE_INDEX(j, k) (k)*(local_jmax) + (j)
#define J_PLANE_INDEX(i, k) (k)*(local_imax) + (i)
#define K_PLANE_INDEX(i, j) (j)*(local_imax) + (i)

void SweepDD(int d, Grid_Data *grid_data, double *volume,
             double *sigt, double *sors, double *psi,
             double *i_plane_psi, double *j_plane_psi, double *k_plane_psi,
             double *psi_lf, double *psi_fr, double *psi_bo)
{
  int i, j, k, local_imax, local_jmax, local_kmax;
  int local_imax_1, local_jmax_1;
  Directions *directions = grid_data->directions;
  int id, jd, kd, istartz, jstartz, kstartz;
  int istopz, jstopz, kstopz;
  int in, jn, kn, im, jm, km;
  int il, ir, jf, jb, kb, kt;

  double dxnow, dynow, dznow;
  double xcos, ycos, zcos, *dx, *dy, *dz;
  double xcos_dxi, ycos_dyj, zcos_dzk;
  double dxi, dyj, dzk;
  double TWO=2.0;
  double psi_new;

  id = directions[d].id;
  jd = directions[d].jd;
  kd = directions[d].kd;
  xcos = directions[d].xcos;
  ycos = directions[d].ycos;
  zcos = directions[d].zcos;

  local_imax = grid_data->nzones[0];
  local_jmax = grid_data->nzones[1];
  local_kmax = grid_data->nzones[2];
  local_imax_1 = local_imax + 1;
  local_jmax_1 = local_jmax + 1;

  dx = grid_data->deltas[0];
  dy = grid_data->deltas[1];
  dz = grid_data->deltas[2];

  if(id > 0){
    istartz = 0; istopz = local_imax-1; in = 1; il = 0; ir = 1;
  }
  else {
    istartz = local_imax-1; istopz = 0; in = -1; il = 1; ir = 0;
  }

  if(jd > 0){
    jstartz = 0; jstopz = local_jmax-1; jn = 1; jf = 0; jb = 1;
  }
  else {
    jstartz = local_jmax-1; jstopz = 0; jn = -1; jf = 1; jb = 0;
  }

  if(kd > 0){
    kstartz = 0; kstopz = local_kmax-1; kn =  1; kb = 0; kt = 1;
  }
  else {
    kstartz = local_kmax-1; kstopz = 0; kn = -1; kb = 1; kt = 0;
  }

  /* Copy the angular fluxes incident upon this subdomain */
  for(k=0; k<local_kmax; k++){
    for(j=0; j<local_jmax; j++){
      /* psi_lf has length (local_imax+1)*local_jmax*local_kmax */
      psi_lf[Left_INDEX(istartz+il, j, k)] = i_plane_psi[I_PLANE_INDEX(j, k)];
    }
  }

  for(k=0; k<local_kmax; k++){
    for(i=0; i<local_imax; i++){
      /* psi_fr has length local_imax*(local_jmax+1)*local_kmax */
      psi_fr[Front_INDEX(i, jstartz+jf,
                         k)] = j_plane_psi[J_PLANE_INDEX(i, k)];
    }
  }

  for(j=0; j<local_jmax; j++){
    for(i=0; i<local_imax; i++){
      /* psi_bo has length local_imax*local_jmax*(local_kmax+1) */
      psi_bo[Bottom_INDEX(i, j, kstartz+
                          kb)] = k_plane_psi[K_PLANE_INDEX(i, j)];
    }
  }

  /*  Perform transport sweep of the grid 1 cell at a time.   */

  for(k=kstartz; std::labs(k-kstartz)<local_kmax; k+=kn){
    dzk = dz[k+1];
    zcos_dzk = TWO*zcos/dzk;
    for(j=jstartz; std::labs(j-jstartz)<local_jmax; j+=jn){
      dyj = dy[j+1];
      ycos_dyj = TWO*ycos/dyj;
      for(i=istartz; std::labs(i-istartz)<local_imax; i+=in){
        dxi = dx[i+1];
        xcos_dxi = TWO*xcos/dxi;

        /* Calculate new zonal flux */
        psi_new =
          (sors[Zonal_INDEX(i, j, k)]
           + psi_lf[Left_INDEX(i+il, j, k    )]*xcos_dxi
           + psi_fr[Front_INDEX(i, j+jf, k    )]*ycos_dyj
           + psi_bo[Bottom_INDEX(i, j, k+kb )]*zcos_dzk)/
          (xcos_dxi + ycos_dyj + zcos_dzk + sigt[Zonal_INDEX(i, j, k)]);
        psi[Zonal_INDEX(i, j, k)] = psi_new;
        /* Apply diamond-difference relationships */
        psi_lf[Left_INDEX(i+ir, j, k )] = TWO*psi_new -
                                          psi_lf[Left_INDEX(i+il, j, k )];
        psi_fr[Front_INDEX(i, j+jb, k )] = TWO*psi_new -
                                           psi_fr[Front_INDEX(i, j+jf, k )];
        psi_bo[Bottom_INDEX(i, j, k+kt )] =  TWO*psi_new -
                                            psi_bo[Bottom_INDEX(i, j, k+kb )];
      }
    }
  }

  /* Copy the angular fluxes exiting this subdomain */
  for(k=0; k<local_kmax; k++){
    for(j=0; j<local_jmax; j++){
      i_plane_psi[I_PLANE_INDEX(j, k)] = psi_lf[Left_INDEX(istopz+ir, j, k)];
    }
  }

  for(k=0; k<local_kmax; k++){
    for(i=0; i<local_imax; i++){
      j_plane_psi[J_PLANE_INDEX(i, k)] = psi_fr[Front_INDEX(i, jstopz+jb, k)];
    }
  }

  for(j=0; j<local_jmax; j++){
    for(i=0; i<local_imax; i++){
      k_plane_psi[K_PLANE_INDEX(i,
                                j)] = psi_bo[Bottom_INDEX(i, j, kstopz+kt)];
    }
  }

}
