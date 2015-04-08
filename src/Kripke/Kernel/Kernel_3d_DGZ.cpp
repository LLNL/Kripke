#include<Kripke/Kernel/Kernel_3d_DGZ.h>
#include<Kripke/Grid.h>
#include<Kripke/SubTVec.h>


#define LG_PRINT_INFO__
#define KRIPKE_USE_ZONE_SLICES

#ifdef KRIPKE_USE_ESSL
extern "C"
{
void dgemm_ (const char& trans,   const char& transb,
                         const int& m1,       const int& n,
                         const int& k,        const double& alpha,
                         const double* a,     const int& lda,
                         const double* b,     const int& ldb,
                         const double& beta,  double* c, const int& ldc);
}
#endif



Kernel_3d_DGZ::Kernel_3d_DGZ() {

}

Kernel_3d_DGZ::~Kernel_3d_DGZ() {

}

Nesting_Order Kernel_3d_DGZ::nestingPsi(void) const {
  return NEST_DGZ;
}

Nesting_Order Kernel_3d_DGZ::nestingPhi(void) const {
  return NEST_DGZ;
}

Nesting_Order Kernel_3d_DGZ::nestingSigt(void) const {
  return NEST_DGZ;
}

Nesting_Order Kernel_3d_DGZ::nestingEll(void) const {
  return NEST_ZGD;
}

Nesting_Order Kernel_3d_DGZ::nestingEllPlus(void) const {
  return NEST_ZDG;
}

Nesting_Order Kernel_3d_DGZ::nestingSigs(void) const {
  return NEST_DGZ;
}


void Kernel_3d_DGZ::LTimes(Grid_Data *grid_data) {
  // Outer parameters
  int nidx = grid_data->total_num_moments;


  for(int ds = 0;ds < grid_data->num_zone_sets;++ ds){
    grid_data->phi[ds]->clear(0.0);
  }

  // Loop over Subdomains
  int num_subdomains = grid_data->subdomains.size();
  for (int sdom_id = 0; sdom_id < num_subdomains; ++ sdom_id){
    Subdomain &sdom = grid_data->subdomains[sdom_id];

    // Get dimensioning
    int num_zones = sdom.num_zones;
    int num_local_groups = sdom.num_groups;
    int num_groups = sdom.phi->groups;
    int group0 = sdom.group0;
    int num_local_directions = sdom.num_directions;
    int num_groups_zones = num_local_groups*num_zones;

    /* 3D Cartesian Geometry */
    double *psi_ptr = sdom.psi->ptr();
    double * KRESTRICT ell = sdom.ell->ptr();
    double * KRESTRICT phi = sdom.phi->ptr(group0, 0, 0);
    
#if 1
     double * KRESTRICT psi = psi_ptr;

     #ifdef KRIPKE_USE_ESSL

     double ONE = 1.0;
     double *ell_dgemm = &ell[0];
     double *psi_ptr_dgemm = &psi_ptr[0];
     double *phi_dgemm = &phi[0];
     //double timestart = MPI_Wtime();
     dgemm_ ('N','N',num_groups_zones, nidx, num_local_directions,ONE, psi_ptr_dgemm,num_groups_zones,ell_dgemm, num_directions,ONE, phi_dgemm, num_zones*num_groups);
     //double dgemmtime = MPI_Wtime() - timestart;

//      #ifdef LG_PRINT_INFO
//	printf("DGZ:  dgemm:  M= %d N= %d K= %d,  time: %g\n",num_groups_zones, nidx, num_local_directions,dgemmtime);
//     #endif


     #else

     for (int i=0; i < nidx; ++i)
       for (int j=0; j < num_local_directions; ++j)
         for (int k = 0; k < num_groups_zones; ++k)
           phi[num_zones*num_groups*i + k] += ell[num_directions*i + j] * psi[num_groups_zones*j + k];

     #endif

   
#else

    for(int nm_offset = 0;nm_offset < nidx;++nm_offset){
      double * KRESTRICT psi = psi_ptr;
      for (int d = 0; d < num_local_directions; d++) {
        double ell_nm_d = ell[d];

#ifdef KRIPKE_USE_OPENMP
#pragma omp parallel for schedule(static)
#endif
        for(int gz = 0;gz < num_groups_zones; ++ gz){
          phi[gz] += ell_nm_d * psi[gz];
        }
        psi += num_groups_zones;
      }
      ell += num_local_directions;
      phi += num_groups*num_zones;
    }
#endif    
  } // Subdomain
}

void Kernel_3d_DGZ::LPlusTimes(Grid_Data *grid_data) {
  // Outer parameters
  int nidx = grid_data->total_num_moments;

  // Loop over Subdomains
  int num_subdomains = grid_data->subdomains.size();
  for (int sdom_id = 0; sdom_id < num_subdomains; ++ sdom_id){
    Subdomain &sdom = grid_data->subdomains[sdom_id];

    // Get dimensioning
    int num_zones = sdom.num_zones;
    int num_local_groups = sdom.num_groups;
    int num_groups = sdom.phi_out->groups;
    int group0 = sdom.group0;
    int num_local_directions = sdom.num_directions;
    int num_groups_zones = num_local_groups*num_zones;

    sdom.rhs->clear(0.0);

    /* 3D Cartesian Geometry */
    double *phi_out_ptr = sdom.phi_out->ptr(group0, 0, 0);
    double * KRESTRICT ell_plus = sdom.ell_plus->ptr();
    double * KRESTRICT rhs = sdom.rhs->ptr();

   //assume column-major storage
//  for (i=0; i < num_local_directions; ++i)
//   for (j=0; j < nidx; ++j) 
//     for (k = 0; k < num_groups_zones; ++k)
//       RHS[k][i * num_local_groups*num_zones] += ELL_PLUS[j][i] * PHI[k][j * num_groups * num_zones];


#if 1
   #ifdef KRIPKE_USE_ESSL
        double *phi_out_dgemm = phi_out_ptr;
        double *rhs_dgemm = rhs;
        double *ell_plus_dgemm = ell_plus;
        double ONE = 1.0; 
        //double timestart = MPI_Wtime();
        dgemm_ ('N','N',num_groups_zones, num_local_directions, nidx,ONE, 
             phi_out_dgemm, num_groups*num_zones, 
             ell_plus_dgemm, nidx, ONE, 
             rhs_dgemm, num_local_groups*num_zones);
        //double dgemmtime = MPI_Wtime() - timestart;
        //printf("DGZ+:  dgemm:  M= %d N= %d K= %d,  time: %g\n",num_groups_zones, num_local_directions,nidx,dgemmtime);

          
   #else  
        double * KRESTRICT phi_out = phi_out_ptr;
        for (int i=0; i < num_local_directions; ++i){
          for (int j=0; j < nidx; ++j) 
            for (int k = 0; k < num_groups_zones; ++k)
              rhs[i*num_local_groups*num_zones + k] += ell_plus[i*nidx + j] * phi_out[j*num_groups *num_zones + k];
        } 
   #endif

#else
    for (int d = 0; d < num_local_directions; d++) {
      double * KRESTRICT phi_out = phi_out_ptr;

      for(int nm_offset = 0;nm_offset < nidx;++nm_offset){
        double ell_plus_d_nm = ell_plus[nm_offset];

#ifdef KRIPKE_USE_OPENMP
#pragma omp parallel for schedule(static)
#endif
        for(int gz = 0;gz < num_groups_zones; ++ gz){
          rhs[gz] += ell_plus_d_nm * phi_out[gz];
        }
        phi_out += num_groups * num_zones;
      }
      ell_plus += nidx;
      rhs += num_local_groups*num_zones;
    }
#endif

  } // Subdomains
}

/**
  Compute scattering source term phi_out from flux moments in phi.
  phi_out(gp,z,nm) = sum_g { sigs(g, n, gp) * phi(g,z,nm) }

  we are mapping sigs(g,d,z) to mean:
    g=source group
    d=legendre coeff
    z=destination group
*/
void Kernel_3d_DGZ::scattering(Grid_Data *grid_data){
  // Loop over zoneset subdomains
  for(int zs = 0;zs < grid_data->num_zone_sets;++ zs){
    // get the phi and phi out references
    SubTVec &phi = *grid_data->phi[zs];
    SubTVec &phi_out = *grid_data->phi_out[zs];
    SubTVec &sigs0 = *grid_data->sigs[0];
    SubTVec &sigs1 = *grid_data->sigs[1];
    SubTVec &sigs2 = *grid_data->sigs[2];

    // get material mix information
    int sdom_id = grid_data->zs_to_sdomid[zs];
    Subdomain &sdom = grid_data->subdomains[sdom_id];
    int const * KRESTRICT mixed_to_zones = &sdom.mixed_to_zones[0];
    int const * KRESTRICT mixed_material = &sdom.mixed_material[0];
    double const * KRESTRICT mixed_fraction = &sdom.mixed_fraction[0];

    // Zero out source terms
    phi_out.clear(0.0);

    // grab dimensions
    int num_mixed = sdom.mixed_to_zones.size();
    int num_zones = sdom.num_zones;
    int num_groups = phi.groups;
    int num_moments = grid_data->total_num_moments;
    int const * KRESTRICT moment_to_coeff = &grid_data->moment_to_coeff[0];

    double *phi_nm_g = phi.ptr();
    double *phi_out_nm = phi_out.ptr();
    for(int nm = 0;nm < num_moments;++ nm){
      // map nm to n
      int n = moment_to_coeff[nm];
      double *sigs0_n_g = sigs0.ptr() + n*num_groups*num_groups;
      double *sigs1_n_g = sigs1.ptr() + n*num_groups*num_groups;
      double *sigs2_n_g = sigs2.ptr() + n*num_groups*num_groups;

      for(int g = 0;g < num_groups;++ g){
        double *phi_out_nm_gp = phi_out_nm;

        for(int gp = 0;gp < num_groups;++ gp){
          double sigs_n_g_gp[3] = {sigs0_n_g[gp], sigs1_n_g[gp], sigs2_n_g[gp]};

          for(int mix = 0;mix < num_mixed;++ mix){
            int zone = mixed_to_zones[mix];
            int material = mixed_material[mix];
            double fraction = mixed_fraction[mix];
            double sigs_value = sigs_n_g_gp[material];

            phi_out_nm_gp[zone] += sigs_value * phi_nm_g[zone] * fraction;
          }

          phi_out_nm_gp += num_zones;
        }
        phi_nm_g += num_zones;
        sigs0_n_g += num_groups;
        sigs1_n_g += num_groups;
        sigs2_n_g += num_groups;
      }
      phi_out_nm += num_groups * num_zones;
    }
  }
}


/**
 * Add an isotropic source, with flux of 1, to every zone with Region 1
 * (or material 0).
 *
 * Since it's isotropic, we're just adding this to nm=0.
 */
void Kernel_3d_DGZ::source(Grid_Data *grid_data){
  // Loop over zoneset subdomains
  for(int zs = 0;zs < grid_data->num_zone_sets;++ zs){
    // get the phi and phi out references
    SubTVec &phi_out = *grid_data->phi_out[zs];

    // get material mix information
    int sdom_id = grid_data->zs_to_sdomid[zs];
    Subdomain &sdom = grid_data->subdomains[sdom_id];
    int const * KRESTRICT mixed_to_zones = &sdom.mixed_to_zones[0];
    int const * KRESTRICT mixed_material = &sdom.mixed_material[0];
    double const * KRESTRICT mixed_fraction = &sdom.mixed_fraction[0];

    // grab dimensions
    int num_mixed = sdom.mixed_to_zones.size();
    int num_zones = sdom.num_zones;
    int num_groups = phi_out.groups;
    int num_moments = grid_data->total_num_moments;

    double *phi_out_nm0_g = phi_out.ptr();
    for(int g = 0;g < num_groups;++ g){
      for(int mix = 0;mix < num_mixed;++ mix){
        int zone = mixed_to_zones[mix];
        int material = mixed_material[mix];
        double fraction = mixed_fraction[mix];

        if(material == 0){
          phi_out_nm0_g[zone] += 1.0 * fraction;
        }
      }

      phi_out_nm0_g += num_zones;
    }
  }
}


/* Sweep routine for Diamond-Difference */
/* Macros for offsets with fluxes on cell faces */
#define I_PLANE_INDEX(j, k) (k)*(local_jmax) + (j)
#define J_PLANE_INDEX(i, k) (k)*(local_imax) + (i)
#define K_PLANE_INDEX(i, j) (j)*(local_imax) + (i)
#define Zonal_INDEX(i, j, k) (i) + (local_imax)*(j) \
  + (local_imax)*(local_jmax)*(k)

void Kernel_3d_DGZ::sweep(Subdomain *sdom) {
  int num_directions = sdom->num_directions;
  int num_groups = sdom->num_groups;
  int num_zones = sdom->num_zones;

  Directions *direction = sdom->directions;

  int local_imax = sdom->nzones[0];
  int local_jmax = sdom->nzones[1];
  int local_kmax = sdom->nzones[2];

  double *dx = &sdom->deltas[0][0];
  double *dy = &sdom->deltas[1][0];
  double *dz = &sdom->deltas[2][0];

  // Upwind/Downwind face flux data
  SubTVec &i_plane = *sdom->plane_data[0];
  SubTVec &j_plane = *sdom->plane_data[1];
  SubTVec &k_plane = *sdom->plane_data[2];

  // All directions have same id,jd,kd, since these are all one Direction Set
  // So pull that information out now
  Grid_Sweep_Block const &extent = sdom->sweep_block;

  std::vector<double> xcos_dxi_all(local_imax);
  std::vector<double> ycos_dyj_all(local_jmax);
  std::vector<double> zcos_dzk_all(local_kmax);

  int *ii_jj_kk_z_idx = extent.ii_jj_kk_z_idx;
  int *offset         = extent.offset;
  int Nslices         = extent.Nhyperplanes;

#ifdef KRIPKE_USE_ZONE_SLICES__
      int N = 11;  // need to parametrize 
      int i_inc = extent.inc_i;
      int j_inc = extent.inc_j;
      int k_inc = extent.inc_k;
      int i_min, i_max, j_min, j_max, k_min, k_max;
      int counter = 0;

      if ( i_inc == 1){
        i_min = extent.start_i;
        i_max = extent.end_i-1;
      }
      else{
        i_min = extent.end_i + 1;
        i_max = extent.start_i;
      }
      if ( j_inc == 1){
        j_min = extent.start_j;
        j_max = extent.end_j-1;
      }
      else{
        j_min = extent.end_j + 1;
        j_max = extent.start_j;
      }
      if ( k_inc == 1){
        k_min = extent.start_k;
        k_max = extent.end_k-1;
      }
      else{
        k_min = extent.end_k + 1;
        k_max = extent.start_k;
      }
      int ii_tmp = (1 - i_inc)/2*i_max;
      int jj_tmp = (1 - j_inc)/2*j_max;
      int kk_tmp = (1 - k_inc)/2*k_max;

      int ii_jj_kk_z_idx[num_zones*4];


      for (int C = 0; C <=(3*N); ++C){   //for each C we can touch zone["i","j","k"]  as well as "d" and "group"    in parallel
       for (int i = 0; i <= C; ++i){
         for (int j = 0; j <= C; ++j){
            int k = C - i - j; // surface equation i+j+j=C
            //flip if needed

            int ii = ii_tmp + i*i_inc;
            int jj = jj_tmp + j*j_inc;
            int kk = kk_tmp + k*k_inc;

            if (ii <= i_max && jj <= j_max && kk <= k_max && ii >= i_min && jj >= j_min && kk >= k_min){
              ii_jj_kk_z_idx[counter*4] = ii;
              ii_jj_kk_z_idx[counter*4+1] = jj;
              ii_jj_kk_z_idx[counter*4+2] = kk;
              ii_jj_kk_z_idx[counter*4+3] = Zonal_INDEX(ii, jj, kk);
              counter++;
           }
         }
       }
     }



#endif  


#ifdef KRIPKE_USE_OPENMP
#pragma omp parallel for
#endif
  for (int d = 0; d < num_directions; ++d) {
    double xcos = direction[d].xcos;
    double ycos = direction[d].ycos;
    double zcos = direction[d].zcos;

    for (int i = 0; i < local_imax; ++i) {
      double dxi = dx[i + 1];
      xcos_dxi_all[i] = 2.0 * xcos / dxi;
    }

    for (int j = 0; j < local_jmax; ++j) {
      double dyj = dy[j + 1];
      ycos_dyj_all[j] = 2.0 * ycos / dyj;
    }

    for (int k = 0; k < local_kmax; ++k) {
      double dzk = dz[k + 1];
      zcos_dzk_all[k] = 2.0 * zcos / dzk;
    }
  }
#ifdef KRIPKE_USE_OPENMP
#pragma omp parallel for collapse(2)
#endif
  for (int d = 0; d < num_directions; ++d) {
    for (int group = 0; group < num_groups; ++group) {
      double * KRESTRICT psi_d_g = sdom->psi->ptr(group, d, 0);
      double * KRESTRICT rhs_d_g = sdom->rhs->ptr(group, d, 0);
      double * KRESTRICT i_plane_d_g = &i_plane(group, d, 0);
      double * KRESTRICT j_plane_d_g = &j_plane(group, d, 0);
      double * KRESTRICT k_plane_d_g = &k_plane(group, d, 0);
      double * KRESTRICT sigt_g = sdom->sigt->ptr(group, 0, 0);

      //printf("ijk_inc= [%d %d %d]; ijk_start= [%d %d %d]; ijk_end= [%d %d %d];\n",extent.inc_i,extent.inc_j,extent.inc_k,
      // 									extent.start_i,extent.start_j,extent.start_k,
      //									extent.end_i,extent.end_j,extent.end_k);
#ifdef KRIPKE_USE_ZONE_SLICES
//we may loop over slices and in each slice have parallel execution for a number of elements
//will need to synchronize after each slice

      for (int element = 0; element < num_zones; element++){
         int ii    = ii_jj_kk_z_idx[element*4];
         int jj    = ii_jj_kk_z_idx[element*4+1];
         int kk    = ii_jj_kk_z_idx[element*4+2];
         int z_idx = ii_jj_kk_z_idx[element*4+3];

          //printf("C = %d  [ii  jj kk] = %d %d %d z_idx = %d\n",C, ii,jj,kk,z_idx);

          double zcos_dzk = zcos_dzk_all[d][kk];
          double ycos_dyj = ycos_dyj_all[d][jj];
          double xcos_dxi = xcos_dxi_all[d][ii];
          int I_P_I = I_PLANE_INDEX(jj, kk); 
          int J_P_I = J_PLANE_INDEX(ii, kk);
          int K_P_I = K_PLANE_INDEX(ii, jj);

           /* Calculate new zonal flux */
            double psi_d_g_z = (rhs_d_g[z_idx]
                + i_plane_d_g[I_P_I] * xcos_dxi
                + j_plane_d_g[J_P_I] * ycos_dyj
                + k_plane_d_g[K_P_I] * zcos_dzk)
                / (xcos_dxi + ycos_dyj + zcos_dzk  + sigt_g[z_idx]);

            psi_d_g[z_idx] = psi_d_g_z;
            /* Apply diamond-difference relationships */
            i_plane_d_g[I_P_I] = 2.0 * psi_d_g_z  - i_plane_d_g[I_P_I];
            j_plane_d_g[J_P_I] = 2.0 * psi_d_g_z  - j_plane_d_g[J_P_I];
            k_plane_d_g[K_P_I] = 2.0 * psi_d_g_z  - k_plane_d_g[K_P_I];

       }

#else
      for (int k = extent.start_k; k != extent.end_k; k += extent.inc_k) {
        double zcos_dzk = zcos_dzk_all[k];
        for (int j = extent.start_j; j != extent.end_j; j += extent.inc_j) {
          double ycos_dyj = ycos_dyj_all[j];
          int z_idx = Zonal_INDEX(extent.start_i, j, k);
          for (int i = extent.start_i; i != extent.end_i; i += extent.inc_i) {
            double xcos_dxi = xcos_dxi_all[i];

            /* Calculate new zonal flux */
            double psi_d_g_z = (rhs_d_g[z_idx]
                + i_plane_d_g[I_PLANE_INDEX(j, k)] * xcos_dxi
                + j_plane_d_g[J_PLANE_INDEX(i, k)] * ycos_dyj
                + k_plane_d_g[K_PLANE_INDEX(i, j)] * zcos_dzk)
                / (xcos_dxi + ycos_dyj + zcos_dzk
                    + sigt_g[z_idx]);

            psi_d_g[z_idx] = psi_d_g_z;
            /* Apply diamond-difference relationships */
            i_plane_d_g[I_PLANE_INDEX(j, k)] = 2.0 * psi_d_g_z
                - i_plane_d_g[I_PLANE_INDEX(j, k)];
            j_plane_d_g[J_PLANE_INDEX(i, k)] = 2.0 * psi_d_g_z
                - j_plane_d_g[J_PLANE_INDEX(i, k)];
            k_plane_d_g[K_PLANE_INDEX(i, j)] = 2.0 * psi_d_g_z
                - k_plane_d_g[K_PLANE_INDEX(i, j)];


            z_idx += extent.inc_i;
          }
        }
      }
#endif      
    } // group
  } // direction

}


