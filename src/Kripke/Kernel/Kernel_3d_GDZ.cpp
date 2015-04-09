#include<Kripke/Kernel/Kernel_3d_GDZ.h>
#include<Kripke/Grid.h>
#include<Kripke/SubTVec.h>


#define LG_PRINT_INFO__

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


Kernel_3d_GDZ::Kernel_3d_GDZ() {

}

Kernel_3d_GDZ::~Kernel_3d_GDZ() {

}

Nesting_Order Kernel_3d_GDZ::nestingPsi(void) const {
  return NEST_GDZ;
}

Nesting_Order Kernel_3d_GDZ::nestingPhi(void) const {
  return NEST_GDZ;
}

Nesting_Order Kernel_3d_GDZ::nestingSigt(void) const {
  return NEST_DGZ;
}

Nesting_Order Kernel_3d_GDZ::nestingEll(void) const {
  return NEST_ZGD;
}

Nesting_Order Kernel_3d_GDZ::nestingEllPlus(void) const {
  return NEST_ZDG;
}

Nesting_Order Kernel_3d_GDZ::nestingSigs(void) const {
  return NEST_GZD;
}


void Kernel_3d_GDZ::LTimes(Grid_Data *grid_data) {
  // Outer parameters
  int nidx = grid_data->total_num_moments;

  // Clear phi
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
    int group0 = sdom.group0;
    int num_local_directions = sdom.num_directions;

    /* 3D Cartesian Geometry */
    double * KRESTRICT phi = sdom.phi->ptr(group0, 0, 0);
    double * KRESTRICT psi_ptr = sdom.psi->ptr();

#ifdef LG_PRINT_INFO
      printf("Kernel_3d_GDZ::LTimes[%d]: NGS= %d  NDS= %d NlocG= %d nidx= %d Nd= %d Nz= %d Ndir= %d\n",Kernel_3d_GDZ_LTimes_counter,num_group_sets, num_dir_sets, num_local_groups, nidx, num_local_directions, num_zones, num_directions);
      printf("Kernel_3d_GDZ::LTimes[%d]: PHI[%d][%d] = ELL[%d][%d] * PSI[%d][%d]\n",Kernel_3d_GDZ_LTimes_counter, num_zones,nidx, num_local_directions,nidx,  num_zones, num_local_directions);
#endif 

    for (int group = 0; group < num_local_groups; ++group) {
      double * KRESTRICT ell_nm = sdom.ell->ptr();

// assume column-wise storage
// for (i=0; i < nidx; ++i)
//   for (j=0; j < num_local_directions; ++j) 
//     for (k = 0; k < num_zones; ++k)
//       PHI[k][i] = ELL[j][i] * PSI[k][j];

//test
#if 1
        double * test_phi_ptr = phi;
        double test_sum = 0.0;

        double * KRESTRICT psi = psi_ptr;


  #ifdef KRIPKE_USE_ESSL
        double ONE = 1.0;
        double *ell_dgemm = &ell_nm[0];
        double *psi_ptr_dgemm = &psi_ptr[0];
        double *phi_dgemm = &phi[0];
        dgemm_ ('N','N',num_zones, nidx, num_local_directions,ONE, psi_ptr_dgemm,num_zones,ell_dgemm, num_directions,ONE, phi_dgemm, num_zones);
  #else
        for (int i=0; i < nidx; ++i)
          for (int j=0; j < num_local_directions; ++j){ 
#ifdef KRIPKE_USE_OPENMP
#pragma omp parallel for
#endif
            for (int k = 0; k < num_zones; ++k)
              phi[num_zones*i + k] += ell_nm[num_local_directions*i + j] * psi_ptr[num_zones*j + k];

          }

  #endif

        
        psi_ptr += num_zones*num_local_directions;
        phi += num_zones*nidx; 

#else
      double * test_phi_ptr = phi;   
      double test_sum = 0.0;
      for(int nm_offset = 0;nm_offset < nidx;++nm_offset){
        double * KRESTRICT psi = psi_ptr;

        for (int d = 0; d < num_local_directions; d++) {
          double ell_nm_d = ell_nm[d];

#ifdef KRIPKE_USE_OPENMP
#pragma omp parallel for
#endif
          for(int z = 0;z < num_zones; ++ z){
            phi[z] += ell_nm_d * psi[z];
          }

          psi += num_zones;
        }
        ell_nm += num_local_directions;
        phi += num_zones;
      }
//test
/*
       for (int i=0; i < nidx; ++i)
            for (int k = 0; k < num_zones; ++k)
              test_sum += test_phi_ptr[num_zones*i + k];

       printf("SUM = %g \n",test_sum);  
*/      
      psi_ptr += num_zones*num_local_directions;
#endif      
    }
  } // Subdomain
}

void Kernel_3d_GDZ::LPlusTimes(Grid_Data *grid_data) {
  // Outer parameters
  int nidx = grid_data->total_num_moments;

  // Loop over Subdomains
  int num_subdomains = grid_data->subdomains.size();
  for (int sdom_id = 0; sdom_id < num_subdomains; ++ sdom_id){
    Subdomain &sdom = grid_data->subdomains[sdom_id];

    // Get dimensioning
    int num_zones = sdom.num_zones;
    int num_local_groups = sdom.num_groups;
    int group0 = sdom.group0;
    int num_local_directions = sdom.num_directions;
#ifdef LG_PRINT_INFO
      printf("Kernel_3d_GDZ::LPlusTimes[%d]: NGS= %d NDS= %d NLG= %d NLD= %d nidx=%d NZ= %d\n",Kernel_3d_GDZ_LPlusTimes_counter, num_group_sets, num_dir_sets, num_local_groups,num_local_directions, nidx, num_zones); 
      printf("Kernel_3d_GDZ::LPlusTimes[%d]: RHS[%d][%d] = ELL[%d][%d] * PHI[%d][%d]\n",Kernel_3d_GDZ_LPlusTimes_counter, num_zones,nidx, nidx,num_local_directions, num_zones, nidx);       
#endif    

    // Get Variables
    sdom.rhs->clear(0.0);

    /* 3D Cartesian Geometry */
    double *ell_plus_ptr = sdom.ell_plus->ptr();

    double * KRESTRICT phi_out_ptr = sdom.phi_out->ptr(group0, 0, 0);
    double * KRESTRICT rhs = sdom.rhs->ptr();

    for (int group = 0; group < num_local_groups; ++group) {
      double *ell_plus = ell_plus_ptr;

// for (i=0; i < num_local_directions; ++i)
//   for (j=0; j < nidx; ++j) 
//     for (k = 0; k < num_zones; ++k)
//       RHS[k][i] += ELL[j][i] * PHI[k][j];
#if 1

   #ifdef KRIPKE_USE_ESSL
        double *phi_out_dgemm = phi_out_ptr;
        double *rhs_dgemm = rhs;
        double *ell_plus_dgemm = ell_plus;
        double ONE = 1.0; 
        dgemm_ ('N','N',num_zones, num_local_directions, nidx,ONE, phi_out_dgemm, num_zones, ell_plus_dgemm, nidx,ONE, rhs_dgemm, num_zones);
   #else
        double * KRESTRICT phi_out = phi_out_ptr;
        for (int i=0; i < num_local_directions; ++i){
          for (int j=0; j < nidx; ++j) 
            for (int k = 0; k < num_zones; ++k)
              rhs[i*num_zones + k] += ell_plus[i*nidx + j] * phi_out[j*num_zones + k];
        } 
   #endif
        phi_out_ptr += num_zones*nidx;
        rhs += num_zones*num_local_directions;
        ell_plus += nidx*num_local_directions;
#else

      for (int d = 0; d < num_local_directions; d++) {
        double * KRESTRICT phi_out = phi_out_ptr;

        for(int nm_offset = 0;nm_offset < nidx;++nm_offset){
          double ell_plus_d_nm = ell_plus[nm_offset];

#ifdef KRIPKE_USE_OPENMP
#pragma omp parallel for
#endif
          for(int z = 0;z < num_zones; ++ z){
            rhs[z] += ell_plus_d_nm * phi_out[z];
          }
          phi_out += num_zones;
        }
        ell_plus += nidx;
        rhs += num_zones;
      }
      phi_out_ptr += num_zones*nidx;
#endif      
    }
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
void Kernel_3d_GDZ::scattering(Grid_Data *grid_data){
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
    int num_coeff = grid_data->legendre_order+1;
    int const * KRESTRICT moment_to_coeff = &grid_data->moment_to_coeff[0];

    double *phi_g = phi.ptr();
    double *sigs0_g_gp = sigs0.ptr();
    double *sigs1_g_gp = sigs1.ptr();
    double *sigs2_g_gp = sigs2.ptr();
    for(int g = 0;g < num_groups;++ g){

      double *phi_out_gp_nm = phi_out.ptr();
      for(int gp = 0;gp < num_groups;++ gp){

        double *phi_g_nm = phi_g;
        for(int nm = 0;nm < num_moments;++ nm){
          // map nm to n
          int n = moment_to_coeff[nm];
          double sigs_g_gp_n[3] = {sigs0_g_gp[n], sigs1_g_gp[n], sigs2_g_gp[n]};

          for(int mix = 0;mix < num_mixed;++ mix){
            int zone = mixed_to_zones[mix];
            int material = mixed_material[mix];
            double fraction = mixed_fraction[mix];
            double sigs_value = sigs_g_gp_n[material];

            phi_out_gp_nm[zone] += sigs_value * phi_g_nm[zone] * fraction;
          }

          phi_g_nm += num_zones;
          phi_out_gp_nm += num_zones;
        }
        sigs0_g_gp += num_coeff;
        sigs1_g_gp += num_coeff;
        sigs2_g_gp += num_coeff;
      }
      phi_g +=  num_moments*num_zones;
    }
  }
}


/**
 * Add an isotropic source, with flux of 1, to every zone with Region 1
 * (or material 0).
 *
 * Since it's isotropic, we're just adding this to nm=0.
 */
void Kernel_3d_GDZ::source(Grid_Data *grid_data){
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

    double *phi_out_g_nm0 = phi_out.ptr();
    for(int g = 0;g < num_groups;++ g){
      for(int mix = 0;mix < num_mixed;++ mix){
        int zone = mixed_to_zones[mix];
        int material = mixed_material[mix];
        double fraction = mixed_fraction[mix];

        if(material == 0){
          phi_out_g_nm0[zone] += 1.0 * fraction;
        }
      }
      phi_out_g_nm0 += num_moments * num_zones;
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

void Kernel_3d_GDZ::sweep(Subdomain *sdom) {
  int num_directions = sdom->num_directions;
  int num_groups = sdom->num_groups;
  int num_zones = sdom->num_zones;

  Directions *direction = sdom->directions;

  int local_imax = sdom->nzones[0];
  int local_jmax = sdom->nzones[1];
  int local_kmax = sdom->nzones[2];

  double * dx = &sdom->deltas[0][0];
  double * dy = &sdom->deltas[1][0];
  double * dz = &sdom->deltas[2][0];

  // Upwind/Downwind face flux data
  SubTVec &i_plane = *sdom->plane_data[0];
  SubTVec &j_plane = *sdom->plane_data[1];
  SubTVec &k_plane = *sdom->plane_data[2];

  // All directions have same id,jd,kd, since these are all one Direction Set
  // So pull that information out now
  Grid_Sweep_Block const &extent = sdom->sweep_block;
  
//  printf("Kernel_3d_GDZ::sweep: local_[i j l]max = [%d %d %d] extent_range [i,j,k] = [%d : %d  %d : %d  %d : %d] extent.inc [i,j,k] = [%d %d %d]\n",local_imax, local_jmax, local_kmax, extent.start_i, extent.end_i, extent.start_j, extent.end_j, extent.start_k, extent.end_k, extent.inc_i, extent.inc_j, extent.inc_k);

  double xcos_dxi_all[num_directions][local_imax];
  double ycos_dyj_all[num_directions][local_jmax];
  double zcos_dzk_all[num_directions][local_kmax];
  
#ifdef KRIPKE_USE_OPENMP
#pragma omp parallel for
#endif
  for (int d = 0; d < num_directions; ++d) {
      double two_xcos = 2.0 * direction[d].xcos;
      double two_ycos = 2.0 * direction[d].ycos;
      double two_zcos = 2.0 * direction[d].zcos;

      for (int i = 0; i < local_imax; ++i) 
        xcos_dxi_all[d][i] = two_xcos / dx[i+1];
      for (int j = 0; j < local_jmax; ++j) 
        ycos_dyj_all[d][j] = two_ycos / dy[j+1];
      for (int k = 0; k < local_kmax; ++k) 
        zcos_dzk_all[d][k] = two_zcos / dz[k+1];
   }  

#ifdef KRIPKE_USE_OPENMP
#pragma omp parallel for collapse(2)
#endif
  for (int group = 0; group < num_groups; ++group) {
    for (int d = 0; d < num_directions; ++d) {
    
    
      double * KRESTRICT sigt_g = sdom->sigt->ptr(group, 0, 0);
      double * KRESTRICT psi_g_d = sdom->psi->ptr(group, d, 0);
      double * KRESTRICT rhs_g_d = sdom->rhs->ptr(group, d, 0);
      double * KRESTRICT i_plane_g_d = i_plane.ptr(group, d, 0);
      double * KRESTRICT j_plane_g_d = j_plane.ptr(group, d, 0);
      double * KRESTRICT k_plane_g_d = k_plane.ptr(group, d, 0);


      /*  Perform transport sweep of the grid 1 cell at a time.   */
     for (int k = extent.start_k; k != extent.end_k; k += extent.inc_k) {
        double zcos_dzk = zcos_dzk_all[d][k];
        for (int j = extent.start_j; j != extent.end_j; j += extent.inc_j) {    // 1 IOP
          double ycos_dyj = ycos_dyj_all[d][j];

          int z_idx = Zonal_INDEX(extent.start_i, j, k);                        // 5 IOP 
          int I_P_I = I_PLANE_INDEX(j, k);                                      // 2 IOP 


          for (int i = extent.start_i; i != extent.end_i; i += extent.inc_i) {  // 1 IOP
            double xcos_dxi = xcos_dxi_all[d][i];
            int J_P_I = J_PLANE_INDEX(i, k);                                    // 2 IOP
            int K_P_I = K_PLANE_INDEX(i, j);                                    // 2 IOP 

            
/* Calculate new zonal flux */                                      

            double psi_g_d_z = (rhs_g_d[z_idx]                                  // 
                + i_plane_g_d[I_P_I] * xcos_dxi
                + j_plane_g_d[J_P_I] * ycos_dyj
                + k_plane_g_d[K_P_I] * zcos_dzk)
                / (xcos_dxi + ycos_dyj + zcos_dzk  + sigt_g[z_idx]);

            psi_g_d[z_idx] = psi_g_d_z;

            /* Apply diamond-difference relationships */
            i_plane_g_d[I_P_I] = 2.0 * psi_g_d_z - i_plane_g_d[I_P_I];
            j_plane_g_d[J_P_I] = 2.0 * psi_g_d_z - j_plane_g_d[J_P_I];
            k_plane_g_d[K_P_I] = 2.0 * psi_g_d_z - k_plane_g_d[K_P_I];

            z_idx += extent.inc_i;
          }
        }
      }
    }
  }
  // say what we really did
  sweep_mode = SWEEP_SERIAL;
}

