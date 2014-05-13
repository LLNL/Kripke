#include<Kripke/Kernel/Kernel_3d_DGZ.h>
#include<Kripke/User_Data.h>
#include<Kripke/SubTVec.h>
#include<Kripke/LMat.h>

Kernel_3d_DGZ::Kernel_3d_DGZ(){

}

Kernel_3d_DGZ::~Kernel_3d_DGZ(){

}

Nesting_Order Kernel_3d_DGZ::nestingPsi(void) const{
  return NEST_DGZ;
}

Nesting_Order Kernel_3d_DGZ::nestingPhi(void) const{
  return NEST_DGZ;
}


void Kernel_3d_DGZ::scattering(Grid_Data *grid_data){
  int num_moments = grid_data->num_moments;
  int num_groups = grid_data->phi->groups;
  int num_zones = grid_data->num_zones;

  double ***phi_in = grid_data->phi->data;
  double ***phi_out = grid_data->phi_out->data;

  // Begin loop over scattering moments
  int m0 = 0;
  for(int n=0; n < num_moments; n++){
    int num_m = grid_data->ell->numM(n);

    for(int m=0; m < num_m; m++){

      double **phi_in_nm = phi_in[m0+m];
      double **phi_out_nm = phi_out[m0+m];

         // Loop over destination group
         for(int gp=0; gp < num_groups; gp++){

          // Loop over source group

           for(int g=0; g < num_groups; g++){

          // Evaluate sigs  for this (n,g,gp) triplet
          evalSigmaS(grid_data, n, g, gp);

          // Get variables
          double *sig_s = &grid_data->sig_s[0];

          double * __restrict__ phi_out_nm_g = phi_out_nm[g];
          double * __restrict__ phi_in_nm_g = phi_in_nm[g];

          for(int zone=0; zone<num_zones; zone++){
            phi_out_nm_g[zone] += sig_s[zone]*phi_in_nm_g[zone];
          } // z

        } // g
      } // gp


    } // m

    m0 += num_m;
  } // n
}

void Kernel_3d_DGZ::LTimes(Grid_Data *grid_data){
  // Outer parameters
  double ***phi = grid_data->phi->data;
  double ***ell = grid_data->ell->data;
  int num_zones = grid_data->num_zones;
  int num_moments = grid_data->num_moments;

  grid_data->phi->clear(0.0);

  // Loop over Group Sets
  int num_group_sets = grid_data->gd_sets.size();
  for(int gset = 0;gset < num_group_sets;++ gset){
    std::vector<Group_Dir_Set> &dir_sets = grid_data->gd_sets[gset];
    int num_dir_sets = dir_sets.size();

    // Loop over Direction Sets
    for(int dset = 0;dset < num_dir_sets;++ dset){
      Group_Dir_Set &gd_set = dir_sets[dset];

      // Get dimensioning
      int num_local_groups = gd_set.num_groups;
      int group0 = gd_set.group0;
      int num_local_directions = gd_set.num_directions;
      int dir0 = gd_set.direction0;

      // Get Variables
      double ***psi = gd_set.psi->data;

      /* 3D Cartesian Geometry */
      for(int n = 0; n < num_moments; n++){
        double ***phi_n = phi + n*n;
        double **ell_n = ell[n];

        for(int m = -n; m <= n; m++){
          double **phi_nm = phi_n[m+n];
          double *ell_n_m = ell_n[m+n];
          for(int d = 0; d < num_local_directions; d++){
            double **psi_d = psi[d];
            double ell_n_m_d = ell_n_m[d+dir0];

            for(int group = 0; group < num_local_groups; ++group){
              double *  psi_d_g = psi_d[group];
              double *  phi_nm_g = phi_nm[group+group0];
              for(int z = 0; z < num_zones; z++){
                double psi_d_g_z = psi_d_g[z];
                phi_nm_g[z] += ell_n_m_d * psi_d_g_z;
              }
            }
          }
        }
      }

    } // Direction Set
  } // Group Set
}


void Kernel_3d_DGZ::LPlusTimes(Grid_Data *grid_data){
  // Outer parameters
  double ***phi_out = grid_data->phi_out->data;
  double ***ell_plus = grid_data->ell_plus->data;
  int num_zones = grid_data->num_zones;
  int num_moments = grid_data->num_moments;

  // Loop over Group Sets
  int num_group_sets = grid_data->gd_sets.size();
  for(int gset = 0;gset < num_group_sets;++ gset){
    std::vector<Group_Dir_Set> &dir_sets = grid_data->gd_sets[gset];
    int num_dir_sets = dir_sets.size();

    // Loop over Direction Sets
    for(int dset = 0;dset < num_dir_sets;++ dset){
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
      for(int d = 0; d < num_local_directions; d++){
        double **psi_d = rhs[d];
        double **ell_plus_d = ell_plus[d+dir0];

        for(int n = 0; n < num_moments; n++){
          double ***phi_out_n = phi_out + n*n;
          double *ell_plus_d_n = ell_plus_d[n];

          for(int m = 0; m <= 2*n; m++){
            double **phi_out_nm = phi_out_n[m];
            double ell_plus_d_n_m = ell_plus_d_n[m];

            for(int group = 0; group < num_local_groups; ++group){
              double * __restrict__ psi_d_g = psi_d[group];
              double * __restrict__ phi_out_nm_g = phi_out_nm[group+group0];

              for(int z = 0; z < num_zones; z++){
                psi_d_g[z] += ell_plus_d_n_m * phi_out_nm_g[z];
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

void Kernel_3d_DGZ::sweep(Grid_Data *grid_data, Group_Dir_Set *gd_set, double *i_plane_ptr, double *j_plane_ptr, double *k_plane_ptr){
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
  SubTVec i_plane_v(nestingPsi(), num_groups, num_directions, local_jmax*local_kmax, i_plane_ptr);
  SubTVec j_plane_v(nestingPsi(), num_groups, num_directions, local_imax*local_kmax, j_plane_ptr);
  SubTVec k_plane_v(nestingPsi(), num_groups, num_directions, local_imax*local_jmax, k_plane_ptr);
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

  for(int d = 0; d < num_directions; ++d){
    double **psi_d = psi[d];
    double **rhs_d = rhs[d];
    double **psi_lf_d = psi_lf.data[d];
    double **psi_fr_d = psi_fr.data[d];
    double **psi_bo_d = psi_bo.data[d];
    double **psi_internal_all_d = psi_internal_all[d];
    double **i_plane_d = i_plane[d];
    double **j_plane_d = j_plane[d];
    double **k_plane_d = k_plane[d];

    for(int group = 0; group < num_groups; ++group){
      double * __restrict__ psi_d_g = psi_d[group];
      double * __restrict__ rhs_d_g = rhs_d[group];
      double * __restrict__ psi_lf_d_g = psi_lf_d[group];
      double * __restrict__ psi_fr_d_g = psi_fr_d[group];
      double * __restrict__ psi_bo_d_g = psi_bo_d[group];
      double * __restrict__ psi_internal_all_d_g = psi_internal_all_d[group];
      double * __restrict__ i_plane_d_g = i_plane_d[group];
      double * __restrict__ j_plane_d_g = j_plane_d[group];
      double * __restrict__ k_plane_d_g = k_plane_d[group];
      double * __restrict__ sigt_g = sigt[group];


      double xcos = direction[d].xcos;
      double ycos = direction[d].ycos;
      double zcos = direction[d].zcos;

      double *psi_int_lf = psi_internal_all_d_g;
      double *psi_int_fr = psi_internal_all_d_g;
      double *psi_int_bo = psi_internal_all_d_g;


      /* Copy the angular fluxes incident upon this subdomain */
      for(int k=0; k<local_kmax; k++){
        for(int j=0; j<local_jmax; j++){
          /* psi_lf has length (local_imax+1)*local_jmax*local_kmax */
          psi_lf_d_g[Left_INDEX(istartz+il, j,
                                k)] = i_plane_d_g[I_PLANE_INDEX(j, k)];
        }
      }

      for(int k=0; k<local_kmax; k++){
        for(int i=0; i<local_imax; i++){
          /* psi_fr has length local_imax*(local_jmax+1)*local_kmax */
          psi_fr_d_g[Front_INDEX(i, jstartz+jf,
                                 k)] = j_plane_d_g[J_PLANE_INDEX(i, k)];
        }
      }

      for(int j=0; j<local_jmax; j++){
        for(int i=0; i<local_imax; i++){
          /* psi_bo has length local_imax*local_jmax*(local_kmax+1) */
          psi_bo_d_g[Bottom_INDEX(i, j, kstartz+
                                  kb)] = k_plane_d_g[K_PLANE_INDEX(i, j)];
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
            psi_lf_d_g[Left_INDEX(i+il, j, k)]
              += psi_int_lf[Zonal_INDEX(i, j, k)];
            psi_fr_d_g[Front_INDEX(i, j+jf, k)]
              += psi_int_fr[Zonal_INDEX(i, j, k)];
            psi_bo_d_g[Bottom_INDEX(i, j, k+kb)]
              += psi_int_bo[Zonal_INDEX(i, j, k)];

            /* Calculate new zonal flux */
            double psi_d_g_z =
              (rhs_d_g[Zonal_INDEX(i, j, k)]
               + psi_lf_d_g[Left_INDEX(i+il, j, k    )]*xcos_dxi
               + psi_fr_d_g[Front_INDEX(i, j+jf, k    )]*ycos_dyj
               + psi_bo_d_g[Bottom_INDEX(i, j, k+kb )]*zcos_dzk)/
              (xcos_dxi + ycos_dyj + zcos_dzk + sigt_g[Zonal_INDEX(i, j, k)]);
            psi_d_g[Zonal_INDEX(i, j, k)] = psi_d_g_z;
            /* Apply diamond-difference relationships */
            psi_lf_d_g[Left_INDEX(i+ir, j, k )] =
              2.0*psi_d_g_z - psi_lf_d_g[Left_INDEX(i+il, j, k)];
            psi_fr_d_g[Front_INDEX(i, j+jb, k )] =
              2.0*psi_d_g_z - psi_fr_d_g[Front_INDEX(i, j+jf, k )];
            psi_bo_d_g[Bottom_INDEX(i, j, k+kt )] =
              2.0*psi_d_g_z - psi_bo_d_g[Bottom_INDEX(i, j, k+kb)];
          }
        }
      }

      /* Copy the angular fluxes exiting this subdomain */
      for(int k=0; k<local_kmax; k++){
        for(int j=0; j<local_jmax; j++){
          i_plane_d_g[I_PLANE_INDEX(j, k)] =
            psi_lf_d_g[Left_INDEX(istopz+ir, j, k)];
        }
      }

      for(int k=0; k<local_kmax; k++){
        for(int i=0; i<local_imax; i++){
          j_plane_d_g[J_PLANE_INDEX(i, k)] =
            psi_fr_d_g[Front_INDEX(i, jstopz+jb, k)];
        }
      }

      for(int j=0; j<local_jmax; j++){
        for(int i=0; i<local_imax; i++){
          k_plane_d_g[K_PLANE_INDEX(i, j)] =
            psi_bo_d_g[Bottom_INDEX(i, j, kstopz+kt)];
        }
      }

    } // group
  } // direction

}

