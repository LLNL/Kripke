#include<Kripke/Kernel/Kernel_3d_ZDG.h>
#include<Kripke/Grid.h>
#include<Kripke/SubTVec.h>

#include <sys/time.h>  
#include <stdio.h> 

#include "Kripke/cu_utils.h"
#ifdef KRIPKE_USE_OPENMP
#include <omp.h>
#endif


#ifdef KRIPKE_USE_HPCLIB
#include <libhpc.h>
#endif

//#define LPlusTimes_sweep_combined

//#define CPU_TIMER


#ifdef KRIPKE_USE_CUDA
int cuda_LTimes_ZDG(double *d_phi, double *h_psi, double *d_ell,
                    int num_zones, int num_groups, int num_local_directions, int num_local_groups, int nidx);

int cuda_LPlusTimes_ZDG(double *rhs, double *phi_out, double *ell_plus,
                    int num_zones, int num_groups, int num_local_directions, int num_local_groups, int nidx);

int cuda_sweep_ZDG( double *rhs, double *phi,
                    double *psi, double *sigt,  Directions *direction,
                    double *i_plane, double *j_plane, double *k_plane,
                    int *ii_jj_kk_z_idx, int *h_offset, int *d_offset, double *dx, double *dy, double *dz,
                    int num_zones, int num_directions, int num_groups,
                    int local_imax, int local_jmax, int local_kmax, int Nslices);

int cuda_LPlusTimes_sweep_ZDG( double *phi_out, double *ell_plus,
                    double *psi, double *sigt,  Directions *direction,
                    double *i_plane, double *j_plane, double *k_plane,
                    int *ii_jj_kk_z_idx, int *h_offset, int *d_offset, double *dx, double *dy, double *dz,
                    int num_zones, int num_directions, int num_groups,
                    int local_imax, int local_jmax, int local_kmax, int Nslices, int nidx);
										
					
int  cuda_scattering_ZDG(int *d_mixed_to_zones, int *d_mixed_material, double *d_mixed_fraction, int *d_mixed_offset,
                      double *d_phi, double *d_phi_out, double *d_sigs0, double *d_sigs1, double *d_sigs2, 
                      int *moment_to_coeff, int num_mixed, int num_moments, int num_groups);

int  cuda_source_ZDG(int *d_mixed_to_zones, int *d_mixed_material, double *d_mixed_fraction, int *d_mixed_offset,
                      double *d_phi_out, int num_mixed, int num_moments, int num_groups);

#endif





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

Nesting_Order Kernel_3d_ZDG::nestingSigt(void) const {
  return NEST_DZG;
}

Nesting_Order Kernel_3d_ZDG::nestingEll(void) const {
  return NEST_ZDG;
}

Nesting_Order Kernel_3d_ZDG::nestingEllPlus(void) const {
  return NEST_ZDG;
}

Nesting_Order Kernel_3d_ZDG::nestingSigs(void) const {
  return NEST_DGZ;
}


void Kernel_3d_ZDG::LTimes(Grid_Data *grid_data) {
  // Outer parameters
  int nidx = grid_data->total_num_moments;

  // Clear phi
  for(int ds = 0;ds < grid_data->num_zone_sets;++ ds){
    grid_data->phi[ds]->clear(0.0);
#ifdef KRIPKE_USE_CUDA
    set_cudaMemZeroAsync( (void *) grid_data->d_phi[ds], (size_t)(grid_data->phi[ds]->elements) * sizeof(double));
#else
    grid_data->phi[ds]->clear(0.0);
#endif

  }

  int num_subdomains = grid_data->subdomains.size();

  // Loop over Subdomains
  for (int sdom_id = 0; sdom_id < num_subdomains; ++ sdom_id){
    Subdomain &sdom = grid_data->subdomains[sdom_id];

    // Get dimensioning
    int num_groups = sdom.phi->groups;
    int num_zones = sdom.num_zones;
    int num_local_groups = sdom.num_groups;
    int group0 = sdom.group0;
    int num_local_directions = sdom.num_directions;
    

    /* 3D Cartesian Geometry */
    double * KRESTRICT ell_d_ptr = sdom.ell->ptr();

#ifdef KRIPKE_USE_CUDA
  if(sweep_mode == SWEEP_GPU){
      cuda_LTimes_ZDG(sdom.d_phi, sdom.psi->ptr(0, 0, 0), sdom.d_ell,
                    num_zones, num_groups, num_local_directions, num_local_groups, nidx);
  }
#endif
 


  if(sweep_mode != SWEEP_GPU){
  int d_unrl = 1;


#ifdef KRIPKE_USE_HPCLIB
   hpmStart(1,"LTimes_ZDG");
#endif

  #ifdef KRIPKE_USE_OPENMP
  #pragma omp parallel 
  { 
     int d; 
  //   int tid = omp_get_thread_num();
  //   double tst, tend;
  //   if (tid==0) tst = omp_get_wtime();
    #pragma omp for
  #endif
#if 0
    for (int z = 0; z < num_zones; z++) {
      double * KRESTRICT psi = sdom.psi->ptr(0, 0, z);
      double * KRESTRICT ell_d = ell_d_ptr;

      for (d = 0; d < num_local_directions-(d_unrl-1); d+=d_unrl) {
        double * KRESTRICT phi = sdom.phi->ptr(group0, 0, z);

        for(int nm_offset = 0;nm_offset < nidx;++nm_offset){
          //double ell_d_nm = ell_d[nm_offset];

          for (int group = 0; group < num_local_groups; ++group) {
            double sum = 0.0;
            for (int ii=0; ii < d_unrl; ++ii) sum += ell_d[nm_offset+(d+ii)*nidx] * psi[group+(d+ii)*num_local_groups];
            phi[group+nm_offset*num_groups] += sum;//    ell_d[nm_offset+d*nidx] * psi[group+d*num_local_groups];
          }
        }
      }
      int dd = d;
      for (int d = dd; d < num_local_directions; ++d) {
        double * KRESTRICT phi = sdom.phi->ptr(group0, 0, z);

        for(int nm_offset = 0;nm_offset < nidx;++nm_offset){
          //double ell_d_nm = ell_d[nm_offset];

          for (int group = 0; group < num_local_groups; ++group) {
            phi[group+nm_offset*num_groups] +=  ell_d[nm_offset+d*nidx] * psi[group+d*num_local_groups];
          }
        }
      }


    }

#else
    for (int z = 0; z < num_zones; z++) {
      double * KRESTRICT psi = sdom.psi->ptr(0, 0, z);
      double * KRESTRICT ell_d = ell_d_ptr;

      for (int d = 0; d < num_local_directions; d++) {
        double * KRESTRICT phi = sdom.phi->ptr(group0, 0, z);

        for(int nm_offset = 0;nm_offset < nidx;++nm_offset){
          double ell_d_nm = ell_d[nm_offset];

          for (int group = 0; group < num_local_groups; ++group) {
            phi[group] += ell_d_nm * psi[group];
          }
          phi += num_groups;
        }
        ell_d += nidx;
        psi += num_local_groups;
      }
    }
#endif
  //    #ifdef KRIPKE_USE_OPENMP
    //  if (tid==0) {
    //    tend = omp_get_wtime();
    //    printf("CPU LTimes_ZDG: time= %g, size(C)= %d, size(ELL)= %d, size(PSI)= %d \n ",tend-tst, sdom.phi->elements, sdom.ell->elements, sdom.psi->elements);
    //  }
   //   #endif
    }
#ifdef BHPMLIB
   hpmStop(1);
#endif

   }
  } // Subdomain

  //keep data on GPU
#ifdef KRIPKE_USE_CUDA__

  if(sweep_mode == SWEEP_GPU){
      for(int ds = 0;ds < grid_data->num_zone_sets;++ ds){
        double *tphi = grid_data->phi[ds]->ptr();
        do_cudaMemcpyD2H ( (void *) tphi, (void*) grid_data->d_phi[ds], grid_data->phi[ds]->elements * sizeof(double));
      }
    }
#endif
}

void Kernel_3d_ZDG::LPlusTimes(Grid_Data *grid_data) {
  // Outer parameters
#ifdef LPlusTimes_sweep_combined
  return;
#endif

  int nidx = grid_data->total_num_moments;

  // Loop over Subdomains
  int num_subdomains = grid_data->subdomains.size();
  for (int sdom_id = 0; sdom_id < num_subdomains; ++ sdom_id){
    Subdomain &sdom = grid_data->subdomains[sdom_id];

    // Get dimensioning
    int num_zones = sdom.num_zones;
    int num_groups = sdom.phi_out->groups;
    int num_local_groups = sdom.num_groups;
    int group0 = sdom.group0;
    int num_local_directions = sdom.num_directions;
    
    //printf("  LPlusTimes:  gset= %d dset=%d group0= %d dir0= %d\n",gset,dset,group0,dir0);
    // Get Variables

    /* 3D Cartesian Geometry */
    double * KRESTRICT ell_plus_ptr = sdom.ell_plus->ptr();

#ifdef KRIPKE_USE_CUDA
  if(sweep_mode == SWEEP_GPU){
      double *dptr_h_rhs = sdom.rhs->ptr();
      if ( sdom.d_rhs == NULL){ // allocate , consider moving to preprocessing;
         sdom.d_rhs = (double *) get_cudaMalloc((size_t) ( sdom.num_zones * sdom.num_groups * sdom.num_directions) * sizeof(double));
      }

      set_cudaMemZeroAsync( (void *) sdom.d_rhs,  (size_t) ( sdom.num_zones * sdom.num_groups * sdom.num_directions ) * sizeof(double));

      cuda_LPlusTimes_ZDG(sdom.d_rhs,sdom.d_phi_out, sdom.d_ell_plus, 
                    num_zones, num_groups, num_local_directions, num_local_groups, nidx);


       sweep_mode = SWEEP_GPU;
  }
#endif


  if(sweep_mode != SWEEP_GPU){

  sdom.rhs->clear(0.0);

  #ifdef KRIPKE_USE_OPENMP
  #pragma omp parallel for
  #endif
    for (int z = 0; z < num_zones; z++) {
      double * KRESTRICT rhs = sdom.rhs->ptr(0, 0, z);

      double * ell_plus_d = ell_plus_ptr;
      for (int d = 0; d < num_local_directions; d++) {

        double * KRESTRICT phi_out = sdom.phi_out->ptr(group0, 0, z);

        for(int nm_offset = 0;nm_offset < nidx;++nm_offset){
          double ell_plus_d_n_m = ell_plus_d[nm_offset];

          for (int group = 0; group < num_local_groups; ++group) {
            rhs[group] += ell_plus_d_n_m * phi_out[group];
          }
          phi_out += num_groups;
        }
        rhs += num_local_groups;
        ell_plus_d += nidx;
      }
    }
  }//
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
void Kernel_3d_ZDG::scattering(Grid_Data *grid_data){
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
 
    static int ER_COUNTER=0;

    // Zero out source terms

#ifdef KRIPKE_USE_CUDA
 if(sweep_mode == SWEEP_GPU)
     set_cudaMemZeroAsync( (void *) grid_data->d_phi_out[zs], (size_t)(grid_data->phi_out[zs]->elements) * sizeof(double));
  else
     phi_out.clear(0.0);
#else
     phi_out.clear(0.0);
#endif

    // grab dimensions
    int num_mixed = sdom.mixed_to_zones.size();
    int num_zones = sdom.num_zones;
    int num_groups = phi.groups;
    int num_coeff = grid_data->legendre_order+1;
    int num_moments = grid_data->total_num_moments;
    int const * KRESTRICT moment_to_coeff = &grid_data->moment_to_coeff[0];

    double *sigs[3] = {
        sigs0.ptr(),
        sigs1.ptr(),
        sigs2.ptr()
    };
   
  
#ifdef KRIPKE_USE_CUDA
 if(sweep_mode == SWEEP_GPU){
  if (sdom.d_mixed_offset == NULL){

    int Nwarps, max_Nwarps = 480;
    if (num_mixed < max_Nwarps)
      Nwarps = num_mixed; 
    else
      Nwarps = max_Nwarps; 

    int * offsets =  sdom.mixed_offset;
    int chunk_size = (num_mixed + Nwarps - 1)/Nwarps;
    for (int i = 0; i < Nwarps; ++i){
      offsets[i] = i*chunk_size;
      if (offsets[i] > num_mixed)  offsets[i] = num_mixed;
    }
    offsets[Nwarps] = num_mixed;
  
    for (int i = 1; i < Nwarps; ++i){
      int shift_back = 0;
      //zone_start
      int zone = mixed_to_zones[ offsets[i] ];
      int zone_m1 = mixed_to_zones[ offsets[i] - 1];
      while (zone == zone_m1){
        zone_m1 = mixed_to_zones[ offsets[i] - 1 - (shift_back+1) ]; 
        shift_back++;
    //    printf("warp_id = %d, zone = %d, zone_m1 = %d, shift_back = %d\n",i,zone, zone_m1, shift_back);
      }
      offsets[i] = offsets[i] - shift_back;
    } 
    for (int i = Nwarps; i <= max_Nwarps; i++)
       offsets[i+1] = offsets[Nwarps];
  
   // for (int i = 0; i <= Nwarps; ++i)
   //    printf("offsets[%d] = %d\n",i,offsets[i]);

//LG need to add validation, to make sure that blocks do not overlap
    for (int i = 1; i < Nwarps; i++){
      if (offsets[i] < offsets[i-1]) 
        printf("ERROR:  offsets are not properly set up\n" );
    }

    sdom.d_mixed_offset = (int*) get_cudaMalloc(size_t   (max_Nwarps+1) * sizeof(int)   );
    do_cudaMemcpyH2D( (void *) (sdom.d_mixed_offset), (void *) offsets, (size_t) (max_Nwarps+1) * sizeof(int));
  }

     double *d_sigs[3];
     d_sigs[0] = grid_data->d_sigs[0];
     d_sigs[1] = grid_data->d_sigs[1];
     d_sigs[2] = grid_data->d_sigs[2];

     cuda_scattering_ZDG(sdom.d_mixed_to_zones, sdom.d_mixed_material, sdom.d_mixed_fraction, sdom.d_mixed_offset,
                      sdom.d_phi, sdom.d_phi_out, grid_data->d_sigs[0],  grid_data->d_sigs[1], grid_data->d_sigs[2],   grid_data->d_moment_to_coeff,  
                      sdom.mixed_to_zones.size(), grid_data->total_num_moments, phi.groups);
                    
}
#endif

 if(sweep_mode != SWEEP_GPU){ 
 
  #if KRIPKE_USE_OPENMP
  int chunk_size[omp_get_max_threads()];
  int offsets[omp_get_max_threads()];
  #pragma omp parallel
  {
  int nthreads = omp_get_num_threads();
  int thread_id = omp_get_thread_num();
  offsets[thread_id] = 0;
  chunk_size[thread_id] = (num_mixed + nthreads - 1)/nthreads;
   
  #pragma omp barrier
  for (int i = 1; i < nthreads; ++i)
    offsets[thread_id] = offsets[thread_id] + chunk_size[thread_id-1];

  //make sure that processing the same zone is done by the same thread.

   int shift_back = 0;
  //zone_start
  if (thread_id > 0){
    int zone = mixed_to_zones[ offsets[thread_id] ];
    int zone_m1 = mixed_to_zones[ offsets[thread_id] - 1];
    while (zone == zone_m1){
      zone_m1 = mixed_to_zones[ offsets[thread_id] - 1 - (shift_back+1) ]; 
      shift_back++;
    }
  } 
  #pragma omp barrier   
  offsets[thread_id] = offsets[thread_id] - shift_back;
//  printf("scatter: org: offset[%d]=  %d, new offset[%d]= %d\n",thread_id,
//	   offsets[thread_id] + shift_back, thread_id,offsets[thread_id]);
     
//  #pragma omp barrier
  //recompute chunk size
//  if (thread_id < (nthreads-1))
//    chunk_size[thread_id] = offsets[thread_id+1] - offsets[thread_id];    
//  else
//    chunk_size[thread_id] = num_mixed - offsets[thread_id];     
}
#endif


    
  int mix_min = 0;
  int mix_max = num_mixed;
#if KRIPKE_USE_OPENMP
#pragma omp parallel private(mix_min,mix_max)
{
   int thread_id = omp_get_thread_num();
   int nthreads = omp_get_num_threads();

   mix_min = offsets[thread_id];
   if (thread_id <  (nthreads-1))
      mix_max = offsets[thread_id+1];
   else
      mix_max = num_mixed;
#endif

   for(int mix = mix_min;mix < mix_max;++ mix){
      int zone = mixed_to_zones[mix];
      int material = mixed_material[mix];
      double fraction = mixed_fraction[mix];
      double *sigs_mat = sigs[material];
      double *phi_z_nm = phi.ptr() + zone*num_groups*num_moments;
      double *phi_out_z_nm = phi_out.ptr() + zone*num_groups*num_moments;

      for(int nm = 0;nm < num_moments;++ nm){
        // map nm to n
        int n = moment_to_coeff[nm];
        double *sigs_n_g = sigs_mat + n*num_groups*num_groups;

        for(int g = 0;g < num_groups;++ g){
          double *phi_out_z_gp = phi_out.ptr() + zone*num_groups*num_moments; //LG whi needed?

          for(int gp = 0;gp < num_groups;++ gp){
            phi_out_z_nm[gp] += sigs_n_g[gp] * phi_z_nm[g] * fraction;
          }
          sigs_n_g += num_groups;
        }
        phi_z_nm += num_groups;
        phi_out_z_nm += num_groups;
      }
    }
#ifdef KRIPKE_USE_OPENMP
 }
#endif
   }
  }

}

/**
 * Add an isotropic source, with flux of 1, to every zone with Region 1
 * (or material 0).
 *
 * Since it's isotropic, we're just adding this to nm=0.
 */
void Kernel_3d_ZDG::source(Grid_Data *grid_data){
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
    int num_zones = sdom.num_zones; //no need
    int num_groups = phi_out.groups;
    int num_moments = grid_data->total_num_moments;

#ifdef KRIPKE_USE_CUDA
    if(sweep_mode == SWEEP_GPU){

    cuda_source_ZDG(sdom.d_mixed_to_zones, sdom.d_mixed_material, sdom.d_mixed_fraction, 
                    sdom.d_mixed_offset, sdom.d_phi_out, 
                    sdom.mixed_to_zones.size(), grid_data->total_num_moments, 
                    phi_out.groups);

    //double *tphi_out = grid_data->phi_out[zs]->ptr();
    //do_cudaMemcpyD2H ( (void *) tphi_out, (void*) grid_data->d_phi_out[zs], grid_data->phi_out[zs]->elements * sizeof(double));
    }
#endif


    if(sweep_mode != SWEEP_GPU){

      for(int mix = 0;mix < num_mixed;++ mix){
      //  int zone = mixed_to_zones[mix];//no need to load if material != 0
        int material = mixed_material[mix];
      //  double fraction = mixed_fraction[mix]; //no need to load if material != 0

        if(material == 0){
          int zone = mixed_to_zones[mix];
          double fraction = mixed_fraction[mix];

          double *phi_out_z_nm0 = phi_out.ptr() + zone*num_moments*num_groups;
          for(int g = 0;g < num_groups;++ g){
            phi_out_z_nm0[g] += 1.0 * fraction;
          }
        }
      }
    }//end of  "if(sweep_mode != SWEEP_GPU)"
  }
}

/* Sweep routine for Diamond-Difference */
/* Macros for offsets with fluxes on cell faces */
#define I_PLANE_INDEX(j, k) (k)*(local_jmax) + (j)
#define J_PLANE_INDEX(i, k) (k)*(local_imax) + (i)
#define K_PLANE_INDEX(i, j) (j)*(local_imax) + (i)
#define Zonal_INDEX(i, j, k) (i) + (local_imax)*(j) \
  + (local_imax)*(local_jmax)*(k)


void Kernel_3d_ZDG::sweep(Subdomain *sdom) {

#ifdef LPlusTimes_sweep_combined
  LPlusTimes_sweep(sdom);
  return;
#endif

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
  
  int *ii_jj_kk_z_idx = extent.ii_jj_kk_z_idx;
  int *offset         = extent.offset;
  int Nslices         = extent.Nhyperplanes;


#ifdef KRIPKE_USE_CUDA
  if(sweep_mode == SWEEP_GPU){

  //LG allocate deltas on GPU
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


  //LG allocate directions on GPU
     if ( sdom->d_directions == NULL){
         sdom->d_directions = (Directions*) get_cudaMalloc((size_t)  num_directions * sizeof(Directions) );
         do_cudaMemcpyH2D( (void*) sdom->d_directions , (void *)  sdom->directions,  (size_t)   num_directions * sizeof(Directions) );
     }
     if ( sdom->d_sigt == NULL){
        sdom->d_sigt = (double *) get_cudaMalloc((size_t) (num_zones*num_groups) * sizeof(double));
        do_cudaMemcpyH2D( (void*) sdom->d_sigt,  (void *)  sdom->sigt->ptr(), (size_t) (num_zones*num_groups) * sizeof(double));
     }

     cuda_sweep_ZDG( sdom->d_rhs, sdom->phi->ptr(),
                     sdom->psi->ptr(), sdom->d_sigt,  sdom->d_directions,
                     i_plane.ptr(),j_plane.ptr(),k_plane.ptr(),
                     extent.d_ii_jj_kk_z_idx, offset, extent.d_offset,
                     sdom->d_delta_x, sdom->d_delta_y, sdom->d_delta_z,
                     num_zones, num_directions, num_groups,
                     local_imax,local_jmax, local_kmax,
                     Nslices);

     // say what we really did
     sweep_mode = SWEEP_GPU;
     return;
  }
#endif

  if(sweep_mode == SWEEP_HYPERPLANE){


#ifdef KRIPKE_USE_OPENMP
    #pragma omp parallel
    {
#endif
    for (int slice = 0; slice < Nslices; slice++){
  struct timeval tim;
  double t1,t2;

#ifdef CPU_TIMER
#pragma omp master
    {
    gettimeofday(&tim, NULL);  
    t1=tim.tv_sec+(tim.tv_usec/1000000.0); 
    }
#endif

#ifdef KRIPKE_USE_OPENMP
#pragma omp for
#endif
    for (int element = offset[slice]; element < offset[slice+1]; ++element){  //for (element = blockId.x + offset[slice]; element < offset[slice+1]; ++element){ 
      int i    = ii_jj_kk_z_idx[element*4];
      int j    = ii_jj_kk_z_idx[element*4+1];
      int k    = ii_jj_kk_z_idx[element*4+2];
      int z    = ii_jj_kk_z_idx[element*4+3];
      int I_P_I = I_PLANE_INDEX(j, k);
      int J_P_I = J_PLANE_INDEX(i, k);
      int K_P_I = K_PLANE_INDEX(i, j);
      double two_inv_dzk = 2.0/dz[k + 1];
      double two_inv_dyj = 2.0/dy[j + 1];
      double two_inv_dxi = 2.0/dx[i + 1];
      double * KRESTRICT sigt_z = sdom->sigt->ptr(0, 0, z);

      // LG get pointer to data corresponding to d=0;
      // LG assume stride of "num_groups" between directions 
      double * KRESTRICT psi_lf_z_d = i_plane.ptr(0, 0, I_P_I);
      double * KRESTRICT psi_fr_z_d = j_plane.ptr(0, 0, J_P_I);
      double * KRESTRICT psi_bo_z_d = k_plane.ptr(0, 0, K_P_I);

      double * KRESTRICT psi_z_d = sdom->psi->ptr(0, 0, z);
      double * KRESTRICT rhs_z_d = sdom->rhs->ptr(0, 0, z);


      for (int d = 0; d < num_directions; ++d) {    // for (d = threadidx.y; d < num_directions; d += blockDim.y){
          double xcos = direction[d].xcos;
          double ycos = direction[d].ycos;
          double zcos = direction[d].zcos;

          double zcos_dzk = zcos * two_inv_dzk;
          double ycos_dyj = ycos * two_inv_dyj;
          double xcos_dxi = xcos * two_inv_dxi;

          for (int group = 0; group < num_groups; ++group) {  //for (group = threadidx.y; group < num_groups; group += blockDim.x)
            int gd = group + d*num_groups;
            double psi_lf_z_d_group = psi_lf_z_d[gd];
            double psi_fr_z_d_group = psi_fr_z_d[gd];
            double psi_bo_z_d_group = psi_bo_z_d[gd];

            /* Calculate new zonal flux */
            double psi_z_d_g = (rhs_z_d[gd]  
                + psi_lf_z_d_group * xcos_dxi
                + psi_fr_z_d_group * ycos_dyj
                + psi_bo_z_d_group * zcos_dzk)
                / (xcos_dxi + ycos_dyj + zcos_dzk + sigt_z[group]);

            psi_z_d[gd] = psi_z_d_g;
            //psi_z_d[group] = psi_z_d_g;

            /* Apply diamond-difference relationships */
            psi_lf_z_d[gd] = 2.0 * psi_z_d_g - psi_lf_z_d_group;
            psi_fr_z_d[gd] = 2.0 * psi_z_d_g - psi_fr_z_d_group;
            psi_bo_z_d[gd] = 2.0 * psi_z_d_g - psi_bo_z_d_group;
          }
      }
    }

#ifdef CPU_TIMER
    #pragma omp master
    {
      gettimeofday(&tim, NULL);  
      t2=tim.tv_sec+(tim.tv_usec/1000000.0);  
      printf("Sweep CPU: Nzones = %d,   %g seconds \n", offset[slice+1] - offset[slice],t2-t1);  
    }
#endif

  }//end of "for (slice"
  #ifdef KRIPKE_USE_OPENMP
  }
  #endif

    // say what we really did
    sweep_mode = SWEEP_HYPERPLANE;
    return;
  }

  for (int k = extent.start_k; k != extent.end_k; k += extent.inc_k) {
    double dzk = dz[k + 1];
    for (int j = extent.start_j; j != extent.end_j; j += extent.inc_j) {
      double dyj = dy[j + 1];
      for (int i = extent.start_i; i != extent.end_i; i += extent.inc_i) {
        double dxi = dx[i + 1];

        int z = Zonal_INDEX(i, j, k);
        double * KRESTRICT sigt_z = sdom->sigt->ptr(0, 0, z);

#ifdef KRIPKE_USE_OPENMP
#pragma omp parallel for
#endif
        for (int d = 0; d < num_directions; ++d) {
          double xcos = direction[d].xcos;
          double ycos = direction[d].ycos;
          double zcos = direction[d].zcos;

          double zcos_dzk = 2.0 * zcos / dzk;
          double ycos_dyj = 2.0 * ycos / dyj;
          double xcos_dxi = 2.0 * xcos / dxi;

          double * KRESTRICT psi_z_d = sdom->psi->ptr(0, d, z);
          double * KRESTRICT rhs_z_d = sdom->rhs->ptr(0, d, z);

          double * KRESTRICT psi_lf_z_d = i_plane.ptr(0, d, I_PLANE_INDEX(j, k));
          double * KRESTRICT psi_fr_z_d = j_plane.ptr(0, d, J_PLANE_INDEX(i, k));
          double * KRESTRICT psi_bo_z_d = k_plane.ptr(0, d, K_PLANE_INDEX(i, j));

          for (int group = 0; group < num_groups; ++group) {
            /* Calculate new zonal flux */
            double psi_z_d_g = (rhs_z_d[group]
                + psi_lf_z_d[group] * xcos_dxi
                + psi_fr_z_d[group] * ycos_dyj
                + psi_bo_z_d[group] * zcos_dzk)
                / (xcos_dxi + ycos_dyj + zcos_dzk + sigt_z[group]);

            psi_z_d[group] = psi_z_d_g;

            /* Apply diamond-difference relationships */
            psi_z_d_g *= 2.0;
            psi_lf_z_d[group] = psi_z_d_g - psi_lf_z_d[group];
            psi_fr_z_d[group] = psi_z_d_g - psi_fr_z_d[group];
            psi_bo_z_d[group] = psi_z_d_g - psi_bo_z_d[group];
          }
        }
      }
    }
  }

  // say what we really did
  sweep_mode = SWEEP_SERIAL;
}




void Kernel_3d_ZDG::LPlusTimes_sweep(Subdomain *sdom) {
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
  
  int *ii_jj_kk_z_idx = extent.ii_jj_kk_z_idx;
  int *offset         = extent.offset;
  int Nslices         = extent.Nhyperplanes;


#ifdef KRIPKE_USE_CUDA
  if(sweep_mode == SWEEP_GPU){

  //LG allocate deltas on GPU
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


  //LG allocate directions on GPU
     if ( sdom->d_directions == NULL){
         sdom->d_directions = (Directions*) get_cudaMalloc((size_t)  num_directions * sizeof(Directions) );
         do_cudaMemcpyH2D( (void*) sdom->d_directions , (void *)  sdom->directions,  (size_t)   num_directions * sizeof(Directions) );
     }
     if ( sdom->d_sigt == NULL){
        sdom->d_sigt = (double *) get_cudaMalloc((size_t) (num_zones*num_groups) * sizeof(double));
        do_cudaMemcpyH2D( (void*) sdom->d_sigt,  (void *)  sdom->sigt->ptr(), (size_t) (num_zones*num_groups) * sizeof(double));
     }

	 
     cuda_LPlusTimes_sweep_ZDG( sdom->d_phi_out, sdom->d_ell_plus,
                     sdom->psi->ptr(), sdom->d_sigt,  sdom->d_directions,
                     i_plane.ptr(),j_plane.ptr(),k_plane.ptr(),
                     extent.d_ii_jj_kk_z_idx, offset, extent.d_offset,
                     sdom->d_delta_x, sdom->d_delta_y, sdom->d_delta_z,
                     num_zones, num_directions, num_groups,
                     local_imax,local_jmax, local_kmax,
                     Nslices, sdom->ell_plus->groups);

     // say what we really did
     sweep_mode = SWEEP_GPU;
     return;
  }
#endif

  if(sweep_mode == SWEEP_HYPERPLANE){


#ifdef KRIPKE_USE_OPENMP
    #pragma omp parallel
    {
#endif
    for (int slice = 0; slice < Nslices; slice++){
  struct timeval tim;
  double t1,t2;

#ifdef CPU_TIMER
#pragma omp master
    {
    gettimeofday(&tim, NULL);  
    t1=tim.tv_sec+(tim.tv_usec/1000000.0); 
    }
#endif

    double * KRESTRICT ell_plus_ptr = sdom->ell_plus->ptr();

    double rhs_local[num_groups];
    int nidx = sdom->ell_plus->groups;
	  
#ifdef KRIPKE_USE_OPENMP
#pragma omp for
#endif
    for (int element = offset[slice]; element < offset[slice+1]; ++element){  //for (element = blockId.x + offset[slice]; element < offset[slice+1]; ++element){ 
      int i    = ii_jj_kk_z_idx[element*4];
      int j    = ii_jj_kk_z_idx[element*4+1];
      int k    = ii_jj_kk_z_idx[element*4+2];
      int z    = ii_jj_kk_z_idx[element*4+3];
      int I_P_I = I_PLANE_INDEX(j, k);
      int J_P_I = J_PLANE_INDEX(i, k);
      int K_P_I = K_PLANE_INDEX(i, j);
      double two_inv_dzk = 2.0/dz[k + 1];
      double two_inv_dyj = 2.0/dy[j + 1];
      double two_inv_dxi = 2.0/dx[i + 1];
      double * KRESTRICT sigt_z = sdom->sigt->ptr(0, 0, z);

      // LG get pointer to data corresponding to d=0;
      // LG assume stride of "num_groups" between directions 
      double * KRESTRICT psi_lf_z_d = i_plane.ptr(0, 0, I_P_I);
      double * KRESTRICT psi_fr_z_d = j_plane.ptr(0, 0, J_P_I);
      double * KRESTRICT psi_bo_z_d = k_plane.ptr(0, 0, K_P_I);

      double * KRESTRICT psi_z_d = sdom->psi->ptr(0, 0, z);
      double * KRESTRICT rhs_z_d = sdom->rhs->ptr(0, 0, z);

	  
	    
	  
      for (int d = 0; d < num_directions; ++d) {    
          double xcos = direction[d].xcos;
          double ycos = direction[d].ycos;
          double zcos = direction[d].zcos;

          double zcos_dzk = zcos * two_inv_dzk;
          double ycos_dyj = ycos * two_inv_dyj;
          double xcos_dxi = xcos * two_inv_dxi;

		 
//	  double * KRESTRICT phi_out = sdom->phi_out->ptr(group0, 0, z);//LG
          double * KRESTRICT phi_out = sdom->phi_out->ptr(0, 0, z);//LG

          double * ell_plus_d = ell_plus_ptr + nidx*d;

          for (int group = 0; group < num_groups; ++group) 
              rhs_local[group] = 0.0;
		  
          for(int nm_offset = 0;nm_offset < nidx;++nm_offset){
            double ell_plus_d_n_m = ell_plus_d[nm_offset];

            for (int group = 0; group < num_groups; ++group) {
              rhs_local[group] += ell_plus_d_n_m * phi_out[group];
            }
            phi_out += num_groups;
          }
		  
		  
          for (int group = 0; group < num_groups; ++group) {  
            int gd = group + d*num_groups;
            double psi_lf_z_d_group = psi_lf_z_d[gd];
            double psi_fr_z_d_group = psi_fr_z_d[gd];
            double psi_bo_z_d_group = psi_bo_z_d[gd];

            /* Calculate new zonal flux */
            double psi_z_d_g = (rhs_local[group]  
                + psi_lf_z_d_group * xcos_dxi
                + psi_fr_z_d_group * ycos_dyj
                + psi_bo_z_d_group * zcos_dzk)
                / (xcos_dxi + ycos_dyj + zcos_dzk + sigt_z[group]);

            psi_z_d[gd] = psi_z_d_g;
            //psi_z_d[group] = psi_z_d_g;

            /* Apply diamond-difference relationships */
            psi_lf_z_d[gd] = 2.0 * psi_z_d_g - psi_lf_z_d_group;
            psi_fr_z_d[gd] = 2.0 * psi_z_d_g - psi_fr_z_d_group;
            psi_bo_z_d[gd] = 2.0 * psi_z_d_g - psi_bo_z_d_group;
          }
      }
    }

#ifdef CPU_TIMER
    #pragma omp master
    {
      gettimeofday(&tim, NULL);  
      t2=tim.tv_sec+(tim.tv_usec/1000000.0);  
      printf("Sweep CPU: Nzones = %d,   %g seconds \n", offset[slice+1] - offset[slice],t2-t1);  
    }
#endif

  }//end of "for (slice"
  #ifdef KRIPKE_USE_OPENMP
  }
  #endif

    // say what we really did
    sweep_mode = SWEEP_HYPERPLANE;
    return;
  }

  for (int k = extent.start_k; k != extent.end_k; k += extent.inc_k) {
    double dzk = dz[k + 1];
    for (int j = extent.start_j; j != extent.end_j; j += extent.inc_j) {
      double dyj = dy[j + 1];
      for (int i = extent.start_i; i != extent.end_i; i += extent.inc_i) {
        double dxi = dx[i + 1];

        int z = Zonal_INDEX(i, j, k);
        double * KRESTRICT sigt_z = sdom->sigt->ptr(0, 0, z);

#ifdef KRIPKE_USE_OPENMP
#pragma omp parallel for
#endif
        for (int d = 0; d < num_directions; ++d) {
          double xcos = direction[d].xcos;
          double ycos = direction[d].ycos;
          double zcos = direction[d].zcos;

          double zcos_dzk = 2.0 * zcos / dzk;
          double ycos_dyj = 2.0 * ycos / dyj;
          double xcos_dxi = 2.0 * xcos / dxi;

          double * KRESTRICT psi_z_d = sdom->psi->ptr(0, d, z);
          double * KRESTRICT rhs_z_d = sdom->rhs->ptr(0, d, z);

          double * KRESTRICT psi_lf_z_d = i_plane.ptr(0, d, I_PLANE_INDEX(j, k));
          double * KRESTRICT psi_fr_z_d = j_plane.ptr(0, d, J_PLANE_INDEX(i, k));
          double * KRESTRICT psi_bo_z_d = k_plane.ptr(0, d, K_PLANE_INDEX(i, j));

          for (int group = 0; group < num_groups; ++group) {
            /* Calculate new zonal flux */
            double psi_z_d_g = (rhs_z_d[group]
                + psi_lf_z_d[group] * xcos_dxi
                + psi_fr_z_d[group] * ycos_dyj
                + psi_bo_z_d[group] * zcos_dzk)
                / (xcos_dxi + ycos_dyj + zcos_dzk + sigt_z[group]);

            psi_z_d[group] = psi_z_d_g;

            /* Apply diamond-difference relationships */
            psi_z_d_g *= 2.0;
            psi_lf_z_d[group] = psi_z_d_g - psi_lf_z_d[group];
            psi_fr_z_d[group] = psi_z_d_g - psi_fr_z_d[group];
            psi_bo_z_d[group] = psi_z_d_g - psi_bo_z_d[group];
          }
        }
      }
    }
  }

  // say what we really did
  sweep_mode = SWEEP_SERIAL;
}


