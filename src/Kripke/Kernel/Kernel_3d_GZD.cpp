#include<Kripke/Kernel/Kernel_3d_GZD.h>
#include<Kripke/Grid.h>
#include<Kripke/SubTVec.h>

#include "Kripke/cu_utils.h"

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

#ifdef KRIPKE_USE_CUDA


int cuda_LTimes_GZD(double *d_phi, double *h_psi, double *d_ell,
                    int num_groups_zones, int num_local_directions, int nidx, int group0);


int cuda_sweep_GZD( double *d_rhs, double *h_phi, double *h_psi, double *d_sigt, Directions *d_direction,
                    double *h_i_plane, double *h_j_plane, double *h_k_plane,
                    int *d_ii_jj_kk_z_idx, int *h_offset,  int *d_offset, double *d_dx, double *d_dy, double *d_dz,
                    int num_zones, int num_directions, int num_groups,
                    int local_imax, int local_jmax, int local_kmax, int Nslices,
                    int start_i, int start_j, int start_k,
                    int end_i, int end_j, int end_k,
                    int inc_i, int inc_j, int inc_k);
#endif



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

Nesting_Order Kernel_3d_GZD::nestingSigt(void) const {
  return NEST_DGZ;
}

Nesting_Order Kernel_3d_GZD::nestingEll(void) const {
  return NEST_ZDG;
}

Nesting_Order Kernel_3d_GZD::nestingEllPlus(void) const {
  return NEST_ZDG;
}

Nesting_Order Kernel_3d_GZD::nestingSigs(void) const {
  return NEST_GZD;
}


void Kernel_3d_GZD::LTimes(Grid_Data *grid_data) {
  // Outer parameters
  int nidx = grid_data->total_num_moments;

#ifdef KRIPKE_USE_CUDA
//  if(sweep_mode == SWEEP_GPU){
    // Clear phi
    for(int ds = 0;ds < grid_data->num_zone_sets;++ ds)
     set_cudaMemZeroAsync( (void *) grid_data->d_phi[ds], (size_t)(grid_data->phi[ds]->elements) * sizeof(double));
//  }
#else

  // Clear phi
  for(int ds = 0;ds < grid_data->num_zone_sets;++ ds){
    grid_data->phi[ds]->clear(0.0);
  }
#endif

 // Loop over Subdomains
  int num_subdomains = grid_data->subdomains.size();
  for (int sdom_id = 0; sdom_id < num_subdomains; ++ sdom_id){
    Subdomain &sdom = grid_data->subdomains[sdom_id];

    // Get dimensioning
    int num_zones = sdom.num_zones;
    int num_local_groups = sdom.num_groups;
    int num_local_directions = sdom.num_directions;
    int num_groups_zones = num_local_groups*num_zones;
    int group0 = sdom.group0;


#ifdef KRIPKE_USE_CUDA
 
//  if(sweep_mode == SWEEP_GPU){

//LG can be an issue with group0!=0 
//LG need to shift pointer to sdom.d_phi

//      do_cudaMemcpyD2H( (void *) sdom.d_phi,  (void *) sdom.phi->ptr(), (size_t) (sdom.phi->elements * sizeof(double)) );

      cuda_LTimes_GZD(sdom.d_phi, sdom.psi->ptr(), sdom.d_ell,
                    num_groups_zones, num_local_directions, nidx, group0);

      //copy results back to CPU (for now) 
      do_cudaMemcpyD2H( (void *) sdom.phi->ptr(),  (void *) sdom.d_phi, (size_t) (sdom.phi->elements * sizeof(double)) );


//  }
  
#else

  if(sweep_mode != SWEEP_GPU){

     /* 3D Cartesian Geometry */
      double *ell_ptr = sdom.ell->ptr();
      double * KRESTRICT psi = sdom.psi->ptr();
      double * KRESTRICT phi = sdom.phi->ptr(group0, 0, 0);
      double * KRESTRICT ell_d = ell_ptr;



    #ifdef KRIPKE_USE_ESSL
        double ONE = 1.0;
        double *ell_dgemm = &ell_ptr[0];
        double *psi_ptr_dgemm = &psi[0];
        double *phi_dgemm = &phi[0];
        dgemm_ ('N','N', nidx, num_groups_zones, num_local_directions, ONE, ell_dgemm, nidx, psi_ptr_dgemm,num_local_directions, ONE, phi_dgemm, nidx);

    #else

     #if 0
        #ifdef KRIPKE_USE_OPENMP
        #pragma omp parallel for
        #endif
        for (int i = 0; i < num_groups_zones; ++i)
          for (int j = 0; j < num_local_directions; ++j)
            for (int k = 0; k < nidx; ++k)
              phi[i*nidx + k] += ell_d[j*nidx + k] * psi[i*num_local_directions + j];
     #else


        #ifdef KRIPKE_USE_OPENMP
        #pragma omp parallel for
        #endif
        for(int gz = 0;gz < num_groups_zones; ++ gz){
          double * KRESTRICT psi = sdom.psi->ptr() + gz*num_local_directions;
          double * KRESTRICT phi = sdom.phi->ptr(group0, 0, 0);
          double * KRESTRICT ell_d = ell_ptr;

          for (int d = 0; d < num_local_directions; d++) {
            double psi_d = psi[d];

            for(int nm_offset = 0;nm_offset < nidx;++nm_offset){
              phi[nm_offset] += ell_d[nm_offset] * psi_d;
            }
            ell_d += nidx;
          }
        }
    #endif

  #endif // essl
  }
#endif //if def KRIPKE_USE_CUDA
  } // Subdomain


}

void Kernel_3d_GZD::LPlusTimes(Grid_Data *grid_data) {
  // Outer parameters
  int nidx = grid_data->total_num_moments;

  // Loop over Subdomains
  int num_subdomains = grid_data->subdomains.size();
  for (int sdom_id = 0; sdom_id < num_subdomains; ++ sdom_id){
    Subdomain &sdom = grid_data->subdomains[sdom_id];

    // Get dimensioning
    int num_zones = sdom.num_zones;
    int num_local_groups = sdom.num_groups;
    int num_local_directions = sdom.num_directions;
    int num_groups_zones = num_local_groups*num_zones;
    int group0 = sdom.group0;

    sdom.rhs->clear(0.0);

    /* 3D Cartesian Geometry */
    double * KRESTRICT ell_plus_ptr = sdom.ell_plus->ptr();
    
    
#if 1

      double * KRESTRICT rhs = sdom.rhs->ptr(0, 0, 0);
      double * KRESTRICT phi_out = sdom.phi_out->ptr(group0, 0, 0);
      double * KRESTRICT ell_plus_d = ell_plus_ptr;

  #ifdef KRIPKE_USE_ESSL
        double ONE = 1.0;
        double *ell_plus_dgemm = &ell_plus_ptr[0];
        double *rhs_dgemm = &rhs[0];
        double *phi_out_dgemm = &phi_out[0];
        dgemm_ ('T','N',num_local_directions, num_groups_zones, nidx,  ONE, ell_plus_dgemm,nidx, phi_out_dgemm, nidx,  ONE, rhs_dgemm, num_local_directions);
  #else
      for (int i = 0; i < num_groups_zones; ++i)
        for (int j = 0; j < num_local_directions; ++j)
          for (int k = 0; k < nidx; ++k)
            rhs[i*num_local_directions + j] += ell_plus_d[j*nidx + k] * phi_out[i*nidx + k];
  #endif

#else
    

#ifdef KRIPKE_USE_OPENMP
#pragma omp parallel for
#endif
    for(int gz = 0;gz < num_groups_zones; ++ gz){
      double * KRESTRICT rhs = sdom.rhs->ptr() + gz*num_local_directions;
      double * KRESTRICT phi_out = sdom.phi_out->ptr(group0, 0, 0) + gz*nidx;
      double * KRESTRICT ell_plus_d = ell_plus_ptr;

      for (int d = 0; d < num_local_directions; d++) {

        for(int nm_offset = 0;nm_offset < nidx;++nm_offset){
          rhs[d] += ell_plus_d[nm_offset] * phi_out[nm_offset];
        }
        ell_plus_d += nidx;
      }
    }
#endif    
    
  } // Subdomain
}


/**
  Compute scattering source term phi_out from flux moments in phi.
  phi_out(gp,z,nm) = sum_g { sigs(g, n, gp) * phi(g,z,nm) }

  we are mapping sigs(g,d,z) to mean:
    g=source group
    d=legendre coeff
    z=destination group
*/
void Kernel_3d_GZD::scattering(Grid_Data *grid_data){
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
    int num_coeff = grid_data->legendre_order+1;
    int num_moments = grid_data->total_num_moments;
    int const * KRESTRICT moment_to_coeff = &grid_data->moment_to_coeff[0];

    double *phi_g = phi.ptr();
    double *sigs0_g_gp = sigs0.ptr();
    double *sigs1_g_gp = sigs1.ptr();
    double *sigs2_g_gp = sigs2.ptr();
    for(int g = 0;g < num_groups;++ g){

      double *phi_out_gp = phi_out.ptr();
      for(int gp = 0;gp < num_groups;++ gp){

        double *sigs_g_gp[3] = {
          sigs0_g_gp,
          sigs1_g_gp,
          sigs2_g_gp
        };

        for(int mix = 0;mix < num_mixed;++ mix){
          int zone = mixed_to_zones[mix];
          int material = mixed_material[mix];
          double fraction = mixed_fraction[mix];
          double *sigs_g_gp_mat = sigs_g_gp[material];
          double *phi_g_z = phi_g + zone*num_moments;
          double *phi_out_gp_z = phi_out_gp + zone*num_moments;

          for(int nm = 0;nm < num_moments;++ nm){
            // map nm to n
            int n = moment_to_coeff[nm];

            phi_out_gp_z[nm] += sigs_g_gp_mat[n] * phi_g_z[nm] * fraction;
          }
        }
        sigs0_g_gp += num_coeff;
        sigs1_g_gp += num_coeff;
        sigs2_g_gp += num_coeff;
        phi_out_gp += num_zones*num_moments;
      }
      phi_g += num_zones*num_moments;
    }
  }
}


/**
 * Add an isotropic source, with flux of 1, to every zone with Region 1
 * (or material 0).
 *
 * Since it's isotropic, we're just adding this to nm=0.
 */
void Kernel_3d_GZD::source(Grid_Data *grid_data){
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

    double *phi_out_g = phi_out.ptr();
    for(int g = 0;g < num_groups;++ g){
      for(int mix = 0;mix < num_mixed;++ mix){
        int zone = mixed_to_zones[mix];
        int material = mixed_material[mix];
        double fraction = mixed_fraction[mix];

        if(material == 0){
          phi_out_g[zone*num_moments] += 1.0 * fraction;
        }
      }
      phi_out_g += num_moments * num_zones;
    }
  }
}

/* Sweep routine for Diamond-Difference */
/* Macros for offsets with fluxes on cell faces */
#define I_PLANE_INDEX(j, k) ((k)*(local_jmax) + (j))
#define J_PLANE_INDEX(i, k) ((k)*(local_imax) + (i))
#define K_PLANE_INDEX(i, j) ((j)*(local_imax) + (i))
#define Zonal_INDEX(i, j, k) ((i) + (local_imax)*(j) \
  + (local_imax)*(local_jmax)*(k))

void Kernel_3d_GZD::sweep(Subdomain *sdom) {
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


#ifdef KRIPKE_USE_CUDA
  //copy data to d_rhs;

     if (sdom->d_delta_x == NULL){
       sdom->d_delta_x = (double*) get_cudaMalloc(size_t   (local_imax+2) * sizeof(double)   );
       do_cudaMemcpyH2D( (void *) (sdom->d_delta_x), (void *) dx, (size_t) (local_imax+2) * sizeof(double));
     }
     if (sdom->d_delta_y == NULL){
       sdom->d_delta_y = (double*) get_cudaMalloc(size_t   (local_jmax+2) * sizeof(double)   );
       do_cudaMemcpyH2D( (void *) (sdom->d_delta_y), (void *) dy, (size_t) (local_jmax+2) * sizeof(double));
     }
     if (sdom->d_delta_z == NULL){
       sdom->d_delta_z = (double*) get_cudaMalloc(size_t   (local_kmax+2) * sizeof(double)   );
       do_cudaMemcpyH2D( (void *) (sdom->d_delta_z), (void *) dz, (size_t) (local_kmax+2) * sizeof(double));
     }

     if (sdom->two_inv_d_delta_x == NULL){
       sdom->two_inv_d_delta_x = (double*) get_cudaMalloc(size_t   (local_imax+2) * sizeof(double)   );
       double two_inv_d_delta_x[local_imax+2];
       for (int i = extent.start_i; i != extent.end_i; i += extent.inc_i)  two_inv_d_delta_x[i+1] = 2.0 / dx[i + 1];
       do_cudaMemcpyH2D( (void *) (sdom->two_inv_d_delta_x), &two_inv_d_delta_x[0] , (size_t) (local_imax+2) * sizeof(double));
     }
     if (sdom->two_inv_d_delta_y == NULL){
       sdom->two_inv_d_delta_y = (double*) get_cudaMalloc(size_t   (local_jmax+2) * sizeof(double)   );
       double two_inv_d_delta_y[local_jmax+2];
       for (int j = extent.start_j; j != extent.end_j; j += extent.inc_j)  two_inv_d_delta_y[j+1] = 2.0 / dy[j + 1];
       do_cudaMemcpyH2D( (void *) (sdom->two_inv_d_delta_y), &two_inv_d_delta_y[0] , (size_t) (local_jmax+2) * sizeof(double));
     }
     if (sdom->two_inv_d_delta_z == NULL){
       sdom->two_inv_d_delta_z = (double*) get_cudaMalloc(size_t   (local_kmax+2) * sizeof(double)   );
       double two_inv_d_delta_z[local_kmax+2];
       for (int k = extent.start_k; k != extent.end_k; k += extent.inc_k)  two_inv_d_delta_z[k+1] = 2.0 / dz[k + 1];
       do_cudaMemcpyH2D( (void *) (sdom->two_inv_d_delta_z), &two_inv_d_delta_z[0] , (size_t) (local_kmax+2) * sizeof(double));
     }



        
     if ( sdom->d_rhs == NULL){ // allocate , consider moving to preprocessing;
         sdom->d_rhs = (double *) get_cudaMalloc((size_t) ( sdom->num_zones * sdom->num_groups * sdom->num_directions) * sizeof(double));
     }
     do_cudaMemcpyH2D( (void *) sdom->d_rhs, sdom->rhs->ptr(), (size_t) (sdom->num_zones * sdom->num_groups * sdom->num_directions) * sizeof(double)); 

     if ( sdom->d_sigt == NULL){
        sdom->d_sigt = (double *) get_cudaMalloc((size_t) (num_zones*num_groups) * sizeof(double));
        do_cudaMemcpyH2D( (void*) sdom->d_sigt,  (void *)  sdom->sigt->ptr(), (size_t) (num_zones*num_groups) * sizeof(double));
     }
     if ( sdom->d_directions == NULL){
         sdom->d_directions = (Directions*) get_cudaMalloc((size_t)  num_directions * sizeof(Directions) );
         do_cudaMemcpyH2D( (void*) sdom->d_directions , (void *)  sdom->directions,  (size_t)   num_directions * sizeof(Directions) );
     }


     cuda_sweep_GZD(sdom->d_rhs, sdom->phi->ptr(), sdom->psi->ptr(), sdom->d_sigt, 
                    sdom->d_directions,
                    i_plane.ptr(), j_plane.ptr(), k_plane.ptr(),
                    extent.d_ii_jj_kk_z_idx, extent.offset, extent.d_offset, 
                    //sdom->d_delta_x, sdom->d_delta_y, sdom->d_delta_z,
                    sdom->two_inv_d_delta_x, sdom->two_inv_d_delta_y, sdom->two_inv_d_delta_z,
                    num_zones,  num_directions, num_groups,
                    local_imax, local_jmax, local_kmax, extent.Nhyperplanes,
                    extent.start_i, extent.start_j, extent.start_k,
                    extent.end_i, extent.end_j, extent.end_k,
                    extent.inc_i, extent.inc_j, extent.inc_k);

   return;

#endif


  double *i_plane_ptr = i_plane.ptr();
  double *psi_ptr = sdom->psi->ptr();
  double *rhs_ptr = sdom->rhs->ptr();
  double *sigt_ptr = sdom->sigt->ptr();

  int i_plane_zones = local_jmax * local_kmax;
  
 
#ifdef KRIPKE_USE_OPENMP
#pragma omp parallel for
#endif
  for (int group = 0; group < num_groups; ++group) {

//    double *sigt_g = sdom->sigt->ptr(group, 0, 0);
    double *sigt_g = &sigt_ptr[group * num_zones];

    /*  Perform transport sweep of the grid 1 cell at a time.   */
    for (int k = extent.start_k; k != extent.end_k; k += extent.inc_k) {
      double dzk = dz[k + 1];
      double two_dz = 2.0 / dzk;
      for (int j = extent.start_j; j != extent.end_j; j += extent.inc_j) {
        double dyj = dy[j + 1];
        double two_dy = 2.0 / dyj;
        for (int i = extent.start_i; i != extent.end_i; i += extent.inc_i) {
          double dxi = dx[i + 1];
          double two_dx = 2.0 / dxi;

          int z = Zonal_INDEX(i, j, k);

//          double * KRESTRICT psi_g_z = sdom->psi->ptr(group, 0, z);
//          double * KRESTRICT rhs_g_z = sdom->rhs->ptr(group, 0, z);

          double * KRESTRICT psi_g_z = &psi_ptr[group*num_zones*num_directions + z*num_directions];//psi(group, 0, z);
          double * KRESTRICT rhs_g_z = &rhs_ptr[group*num_zones*num_directions + z*num_directions];


          double * KRESTRICT psi_lf_g_z = &i_plane_ptr[group*i_plane_zones*num_directions+I_PLANE_INDEX(j, k)*num_directions];

//          double * KRESTRICT psi_lf_g_z = i_plane.ptr(group, 0, I_PLANE_INDEX(j, k));
          double * KRESTRICT psi_fr_g_z = j_plane.ptr(group, 0, J_PLANE_INDEX(i, k));
          double * KRESTRICT psi_bo_g_z = k_plane.ptr(group, 0, K_PLANE_INDEX(i, j));

          for (int d = 0; d < num_directions; ++d) {
            double xcos = direction[d].xcos;
            double ycos = direction[d].ycos;
            double zcos = direction[d].zcos;

            double zcos_dzk = zcos * two_dz;
            double ycos_dyj = ycos * two_dy;
            double xcos_dxi = xcos * two_dx;

            /* Calculate new zonal flux */
            double psi_g_z_d = (rhs_g_z[d] + psi_lf_g_z[d] * xcos_dxi
                + psi_fr_g_z[d] * ycos_dyj + psi_bo_g_z[d] * zcos_dzk)
                / (xcos_dxi + ycos_dyj + zcos_dzk
                    + sigt_g[Zonal_INDEX(i, j, k)]);

            psi_g_z[d] = psi_g_z_d;

            /* Apply diamond-difference relationships */
            psi_lf_g_z[d] = 2.0 * psi_g_z_d - psi_lf_g_z[d];
            psi_fr_g_z[d] = 2.0 * psi_g_z_d - psi_fr_g_z[d];
            psi_bo_g_z[d] = 2.0 * psi_g_z_d - psi_bo_g_z[d];
          }
        }
      }
    }
  } // group

  // say what we really did
  sweep_mode = SWEEP_SERIAL;
}

void Kernel_3d_GZD::LPlusTimes_sweep(Subdomain *sdom) {
  return ;
}

