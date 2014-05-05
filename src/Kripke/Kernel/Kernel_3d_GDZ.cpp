#include<Kripke/Kernel/Kernel_3d_GDZ.h>
#include<Kripke/user_data.h>
#include<Kripke/SubTVec.h>

Kernel_3d_GDZ::Kernel_3d_GDZ(){

}

Kernel_3d_GDZ::~Kernel_3d_GDZ(){

}

Nesting_Order Kernel_3d_GDZ::nesting(void) const{
  return NEST_GDZ;
}


// Computational Kernels
void Kernel_3d_GDZ::evalSigmaTot(Grid_Data *grid_data){

}
void Kernel_3d_GDZ::evalSigmaS(Grid_Data *grid_data){

}
void Kernel_3d_GDZ::LTimes(Grid_Data *grid_data){

}
void Kernel_3d_GDZ::LPlusTimes(Grid_Data *grid_data){

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

void Kernel_3d_GDZ::sweep(Grid_Data *grid_data, Group_Dir_Set *gd_set){
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

  SubTVec psi_lf(nesting(), num_groups, num_directions,
                  (local_imax+1)*local_jmax*local_kmax);
  SubTVec psi_fr(nesting(), num_groups, num_directions,
                  local_imax*(local_jmax+1)*local_kmax);
  SubTVec psi_bo(nesting(), num_groups, num_directions,
                  local_imax*local_jmax*(local_kmax+1));

  SubTVec i_plane_v(nesting(), num_groups, num_directions, local_jmax*local_kmax);
  SubTVec j_plane_v(nesting(), num_groups, num_directions, local_imax*local_kmax);
  SubTVec k_plane_v(nesting(), num_groups, num_directions, local_imax*local_jmax);
  double ***i_plane = i_plane_v.data;
  double ***j_plane = j_plane_v.data;
  double ***k_plane = k_plane_v.data;

  SubTVec psi_internal(nesting(), num_groups, num_directions, num_zones);
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

  for(int group = 0; group < num_groups; ++group){
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

    for(int d = 0; d < num_directions; ++d){
      double * __restrict__ psi_g_d = psi_g[d];
      double * __restrict__ rhs_g_d = rhs_g[d];
      double * __restrict__ psi_lf_g_d = psi_lf_g[d];
      double * __restrict__ psi_fr_g_d = psi_fr_g[d];
      double * __restrict__ psi_bo_g_d = psi_bo_g[d];
      double * __restrict__ psi_internal_all_g_d = psi_internal_all_g[d];
      double * __restrict__ i_plane_g_d = i_plane_g[d];
      double * __restrict__ j_plane_g_d = j_plane_g[d];
      double * __restrict__ k_plane_g_d = k_plane_g[d];


      double xcos = direction[d].xcos;
      double ycos = direction[d].ycos;
      double zcos = direction[d].zcos;

      double *psi_int_lf = psi_internal_all_g_d;
      double *psi_int_fr = psi_internal_all_g_d;
      double *psi_int_bo = psi_internal_all_g_d;


      /* Copy the angular fluxes incident upon this subdomain */
      for(int k=0; k<local_kmax; k++){
        for(int j=0; j<local_jmax; j++){
          /* psi_lf has length (local_imax+1)*local_jmax*local_kmax */
          psi_lf_g_d[Left_INDEX(istartz+il, j,
                                k)] = i_plane_g_d[I_PLANE_INDEX(j, k)];
        }
      }

      for(int k=0; k<local_kmax; k++){
        for(int i=0; i<local_imax; i++){
          /* psi_fr has length local_imax*(local_jmax+1)*local_kmax */
          psi_fr_g_d[Front_INDEX(i, jstartz+jf,
                                 k)] = j_plane_g_d[J_PLANE_INDEX(i, k)];
        }
      }

      for(int j=0; j<local_jmax; j++){
        for(int i=0; i<local_imax; i++){
          /* psi_bo has length local_imax*local_jmax*(local_kmax+1) */
          psi_bo_g_d[Bottom_INDEX(i, j, kstartz+
                                  kb)] = k_plane_g_d[K_PLANE_INDEX(i, j)];
        }
      }

      /*  Perform transport sweep of the grid 1 cell at a time.   */
      for(int k=kstartz; std::abs(k-kstartz)<local_kmax; k+=kn){
        double dzk = dz[k+1];
        double zcos_dzk = 2.0*zcos/dzk;
        for(int j=jstartz; std::abs(j-jstartz)<local_jmax; j+=jn){
          double dyj = dy[j+1];
          double ycos_dyj = 2.0*ycos/dyj;
          for(int i=istartz; std::abs(i-istartz)<local_imax; i+=in){
            double dxi = dx[i+1];
            double xcos_dxi = 2.0*xcos/dxi;

            /* Add internal surface source data */
            psi_lf_g_d[Left_INDEX(i+il, j,
                                  k)] += psi_int_lf[Zonal_INDEX(i, j, k)];
            psi_fr_g_d[Front_INDEX(i, j+jf,
                                   k)] += psi_int_fr[Zonal_INDEX(i, j, k)];
            psi_bo_g_d[Bottom_INDEX(i, j, k+
                                    kb)] += psi_int_bo[Zonal_INDEX(i, j, k)];

            /* Calculate new zonal flux */
            double psi_g_d_z =
              (rhs_g_d[Zonal_INDEX(i, j, k)]
               + psi_lf_g_d[Left_INDEX(i+il, j, k    )]*xcos_dxi
               + psi_fr_g_d[Front_INDEX(i, j+jf, k    )]*ycos_dyj
               + psi_bo_g_d[Bottom_INDEX(i, j, k+kb )]*zcos_dzk)/
              (xcos_dxi + ycos_dyj + zcos_dzk + sigt_g[Zonal_INDEX(i, j, k)]);
            psi_g_d[Zonal_INDEX(i, j, k)] = psi_g_d_z;
            /* Apply diamond-difference relationships */
            psi_lf_g_d[Left_INDEX(i+ir, j, k )] = 2.0*psi_g_d_z -
                                         psi_lf_g_d[Left_INDEX(i+il,j, k )];
            psi_fr_g_d[Front_INDEX(i, j+jb, k )] = 2.0*psi_g_d_z -
                                         psi_fr_g_d[Front_INDEX(i,j+ jf,k ) ];
            psi_bo_g_d[Bottom_INDEX(i, j, k+kt )] =  2.0*psi_g_d_z -
                                     psi_bo_g_d[Bottom_INDEX(i, j, k + kb )];
          }
        }
      }

      /* Copy the angular fluxes exiting this subdomain */
      for(int k=0; k<local_kmax; k++){
        for(int j=0; j<local_jmax; j++){
          i_plane_g_d[I_PLANE_INDEX(j,k)] =
            psi_lf_g_d[Left_INDEX(istopz+ir, j, k)];
        }
      }

      for(int k=0; k<local_kmax; k++){
        for(int i=0; i<local_imax; i++){
          j_plane_g_d[J_PLANE_INDEX(i, k)] =
            psi_fr_g_d[Front_INDEX(i, jstopz+jb, k)];
        }
      }

      for(int j=0; j<local_jmax; j++){
        for(int i=0; i<local_imax; i++){
          k_plane_g_d[K_PLANE_INDEX(i, j)] =
            psi_bo_g_d[Bottom_INDEX(i, j, kstopz+kt)];
        }
      }

    }
  }

}

