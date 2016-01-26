#include <stdio.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "Directions.h"
#include <Kripke/Subdomain.h>
#include <mpi.h>
#include <omp.h>
#include <Kripke/SubTVec.h>


/*

  DEVELOPMENT OPTIONS - GPU


  Towards fastest DD solution for problem sets which fit in GPU memory.
  ---------------------------------------------------------------------

  0.   Optimize MPI operation such that completed subdomains communicate boundary fluxes immediately.
  
  0.   Measure flux-register kernel bandwidth.
       + fully characterize where any missing bandwidth is coming from.
       + try to include all data transferred from GMEM - and match with nprof measurements
       = 

  1.   Scalability studies.
       + These need to be ongoing as they are the ultimate indicator of overall performance.
       + 1-D studies on surface can currently achieve Manhattan distances of 128
         - larger than the nominal full-scale, 3d, Manhattan distance on Sierra. 
       + can be expanded to 256 just by using the second K40/node
       + can be expanded to 2k (16 cores * 128 nodes ... really 160 nodes?) and using MPS
         - likely config for size would be something like 157 nodes * 14 ranks/node w/ 7ranks/GPU w/ 2 SMs per rank
	   . 72M zones arranged as 416 x 416 x 416 w/ 4 directions and 2 groups. 
	 - likely config for throughput would be something like 108 nodes w/ 2 ranks/node, 1rank/GPU and 15 SMs per rank
	   . 7M zones arranged as 192 x 192 x 192 zones w/ 16 directions and 128 groups 
       + scalabilitie studies needed on both Intel/Nvidia and Power/Nvidia

  1.   Scaling model
       + is the sweep operation still as important/dominant as in the past?
       + what perf/FOM do we expect from Sierra?

  1.   Launch 15 blocks simultaneously - expect one per SMs.  
       + this should minimize launch latency
       + this will also permit using a single stream per subdomain for H2D, Kernels, D2H - help claifiy code.
         - just attach a single stream to each sdom.
       [DONE]

  1.   Cache i-fluxes in SMEM at start/end of kernel.  Currently these fluxes are read from GMEM when needed 
       + this potentially (fairly sure) adds significant latency.

  2.   Cache sigma.  Too big for SMEM?  So probably need some new scheme to make efficient.
       + independent of direction?  So maybe works best when directions used in sets of 16?
       + sigma will likely need to be restrided as zones first?  Need to look at this.  Might be difficult.

  2.   Investigate alternatives for batched GEMM.
       + currently only seeing 122 GF/s using cublasDgemmBatched.  Way off device poential of 1.4 TF/s.
       + Benchmark kernel from CHOLMOD - Darko reports better perf.
       + Investigate benefits of caching L and L+ 
       + can go arbitrarily deep here

  2.   Fold ParticleEdit into sweep.  
       + classic parallel reduction operation
       + Should be essentially free

  3.   General code cleanup.  requires ADAM

  4.   Clean up templatization of octants.
       + current code is easy to debug, but very verbose - makes it look like a lot of CUDA code (which it isn't)

  7.   Generalize handling of large numbers of materials.
       + likely easy if cublasDgemmBatched is successfully replaced

  9.   Generalize the flux register kernel such that it treats groups and directions the same.  
       + this would eliminate the need for directions to be multiples of 4/16
       + this would permit the kernel to be efficiently applied to 'isotropic' problems.
       [NOT REQUIRED?  As soon as the directions are complete, it moves to the next group in the subdomain, 
       which constitutes coalesced access already.]

  10.  Generalize flux register kernel for arbitrary plane dimensions.
       + this would minimize the perf effects when zones are not multiples of 8/16 



  Towards DD of problem sets which don't fit in GPU memory
  --------------------------------------------------------

  0.   Paper studies.  
       + SOL time to copy moment-condensed subdomain data to/from host.  
       + How does this compare to currently measured kernel times?  (Scatter/Source/L+/Sweep/L)
       + potentially investigate opportunities for further data compression
       + gets a bit messy when groups and directions are combined - maybe best to still keep these seperate for now.
       + can only condense one octant worth at a time - so this is only useful when num_moments < num_directions_per_octant?  
         - Unclear.  Need to examine code.
       
  0.   Re-architect Sweep_Solver loop.
       - the GPU operations need to start and stop at (either right before or right after) Scatter/Source
       - this is were the state exists in (typically) most compact form - moments*groups*zones.



  Towards Support for Semi-stuctured 
  ----------------------------------

  0.   Gated by implementation of TriLinear DG in Kripke

  1.   Paper studies / evaluation
       + flop counts vs. bandwidth requirements.  What is the SOL limiter?



  Towards Support for Fully Unstructured
  --------------------------------------

  0.   Gated / completely dependent on support for Semi-structured?

  1.   Paper studies / evaluation


*/



#define KRESTRICT __restrict__

//#define USE_PSI_HOST_MEM

//#define CU_TIMING

#define FluxPlaneDim 8 
#define FluxGroupDim 16
#define NumSMs 15
  
#define MAX ((a<b)?b:a)
#define MIN ((a>b)?b:a)


//Macro for checking cuda errors following a cuda launch or api call
#define cudaCheckError() {                                              \
  cudaError_t e=cudaGetLastError();                                     \
  if(e!=cudaSuccess) {                                                  \
    printf("Cuda failure %s:%d: '%s'\n",__FILE__,__LINE__,cudaGetErrorString(e)); \
    exit(EXIT_FAILURE);                                                 \
  }                                                                     \
}

/* function declarations */

extern "C" {
void dgemm_custom_simple_1block_batch (cudaStream_t stream, 
  		     	               cublasOperation_t transa, 
				       cublasOperation_t transb, 
				       int *mlist, int *nlist, int *klist, 
				       const double *alpha, 
				       const double **Alist, int *ldalist, 
				       const double **Blist, int *ldblist, 
				       const double *beta, 
				       double **Clist, int *ldclist, int nbatch);
}

__global__ void  LTimes_ZGD(double *phi, double * __restrict__ psi, const double * __restrict__ ell,
                                int num_zones, int num_groups, int num_local_directions, int num_local_groups, 
                                int nidx, int group0);

__global__ void  LPlusTimes_ZGD(double *rhs, double * __restrict__ phi_out, double * __restrict__ ell_plus,
                                int num_zones, int num_groups, int num_local_directions, int num_local_groups,
                                int nidx, int group0);

__global__ void scattering_ZGD(int * __restrict__ d_mixed_to_zones, int * __restrict__ d_mixed_material,
                               double * __restrict__ d_mixed_fraction, int * __restrict__ d_mixed_offset,
                               double * __restrict__ d_phi, double *d_phi_out, double * __restrict__ d_sigs0,
                              double * __restrict__ d_sigs1, double * __restrict__ d_sigs2,
                               int * __restrict__ moment_to_coeff, int num_mixed, int num_moments, int num_groups, int num_coeff);

__global__ void scattering_ZGD_step2(const int * __restrict__ mixed_to_zones, const int * __restrict__ mixed_material,
                               const double * __restrict__ mixed_fraction, const int * __restrict__ mixed_offset,
                               const double * __restrict__ phi, double *phi_out, double * __restrict__ sigs0,
                               const double * __restrict__ sigs1, const double * __restrict__ sigs2,
                               const int * __restrict__ moment_to_coeff, int num_mixed, int num_moments, int num_groups, int num_coeff);


__global__ void scattering_ZGD_step3(const int * __restrict__ mixed_to_zones, const int * __restrict__ mixed_material,
                               const double * __restrict__ mixed_fraction, const int * __restrict__ mixed_offset,
                               const double * __restrict__ phi, double *phi_out, double * __restrict__ sigs0,
                               const double * __restrict__ sigs1, const double * __restrict__ sigs2,
                               const int * __restrict__ moment_to_coeff, int num_mixed, int num_moments, int num_groups, int num_coeff);


__global__ void source_ZGD(int * __restrict__ mixed_to_zones, int * __restrict__ mixed_material,
                               double * __restrict__ mixed_fraction, int * __restrict__ mixed_offset,
                               double *phi_out, int num_moments, int num_groups);

__global__ void  sweep_over_hyperplane_ZGD(int sliceID, int * __restrict__ offset, int * __restrict__ ii_jj_kk_z_idx, int num_directions, int num_groups,
                    int local_imax, int local_jmax, int local_kmax, 
                    double * __restrict__ dx, double * __restrict__ dy, double * __restrict__ dz, 
                    double * __restrict__ rhs, double * __restrict__ phi, double * __restrict__ psi, 
                    const double * __restrict__ sigt, const Directions * __restrict__ direction, 
                    double *i_plane, double *j_plane, double *k_plane);


__global__ void LPlusTimes_sweep_over_hyperplane_ZGD(int sliceID, int * __restrict__ offset, int * __restrict__ ii_jj_kk_z_idx, int num_directions, int num_groups, int num_local_groups, int nidx, int group0,
                    int local_imax, int local_jmax, int local_kmax, double * __restrict__ dx, double * __restrict__ dy,
                    double * __restrict__ dz, double *__restrict__ phi_out, double * __restrict__ ell_plus, double * __restrict__ psi,
                    const double * __restrict__ sigt, const Directions * __restrict__ direction,
                    double *i_plane, double *j_plane, double *k_plane);
					
					
					
int cuda_LTimes_ZGD(double *d_phi, double *h_psi, double *d_psi, double *d_ell,
                    int num_zones, int num_groups, int num_local_directions, int num_local_groups, 
                    int nidx, int group0);


int  cuda_LPlusTimes_ZGD(double *d_rhs, double *h_phi_out, double *h_ell_plus,
                    int num_zones, int num_groups, int num_local_directions, int num_local_groups, 
                    int nidx, int group0);

int  cuda_scattering_ZGD(int *d_mixed_to_zones, int *d_mixed_material, double *d_mixed_fraction, int *d_mixed_offset,
                      double *d_phi, double *d_phi_out, double *d_sigs0,
                      double *d_sigs1, double *d_sigs2,
                      int *moment_to_coeff, int num_mixed, int num_moments, int num_groups, int num_coeff);

int  cuda_source_ZGD(int *d_mixed_to_zones, int *d_mixed_material, double *d_mixed_fraction, int *d_mixed_offset, 
                     double *d_phi_out, int num_moments,  int num_groups, int num_mixed);


int cuda_sweep_ZGD( double *rhs, double *phi, double *psi,  double *sigt, Directions *direction,
                    double *i_plane, double *j_plane, double *k_plane,
                    int *ii_jj_kk_z_idx, int *offset, int *d_offset, double *dx, double *dy, double *dz,
                    int num_zones, int num_directions, int num_groups,
                    int local_imax, int local_jmax, int local_kmax, int Nslices);

int cuda_LPlusTimes_sweep_ZGD( double *phi_out, double *ell_plus,
                    double *psi, double *sigt,  Directions *direction,
                    double *i_plane, double *j_plane, double *k_plane,
                    int *ii_jj_kk_z_idx, int *h_offset, int *d_offset, double *dx, double *dy, double *dz,
			       int num_zones, int num_directions, int num_groups, /*int num_local_groups,*/ int nidx,
                    int local_imax, int local_jmax, int local_kmax, int Nslices);


int cuda_scattering_ZGD2 ( Subdomain *sdom, double *d_inflated_materials, int num_moments, int num_groups );


static cublasHandle_t kripke_cbhandle = NULL;

int cuda_scattering_ZGD2 ( Subdomain *sdom, double *d_inflated_material, int num_moments, int num_groups ) 
{

  int num_mixed = sdom->mixed_to_zones.size();

  // cublasBatched routines are much faster for this size of GEMM.
  // Unfortunately they must take a common value for 'alpha' or, in this case, fraction.
  // So, collect all of the GEMMs for zones for which there is only a single material, (alpha=fraction=1) and
  // perform those as a batch.  

  static double ** h_Aarray_scatter = NULL;
  static double ** h_Barray_scatter = NULL;
  static double ** h_Carray_scatter = NULL;
  static double ** d_Aarray_scatter = NULL;
  static double ** d_Barray_scatter = NULL;
  static double ** d_Carray_scatter = NULL;
  int num_zones_scatter = 0;

  if ( ! h_Aarray_scatter ) {

    // batched operation is the same every time - so only need to do this once.
    // small list of pointers to store on device
    // but means that full perf is not achieved until the 3rd iteration
    
    //ref_num_zones = num_zones;
    h_Aarray_scatter = (double **) malloc (3*num_mixed*sizeof(double*));
    h_Barray_scatter = h_Aarray_scatter + num_mixed;
    h_Carray_scatter = h_Barray_scatter + num_mixed;
    
    cudaMalloc ((void**)&d_Aarray_scatter, 3*num_mixed*sizeof(double*) ); 
    d_Barray_scatter = d_Aarray_scatter + num_mixed;
    d_Carray_scatter = d_Barray_scatter + num_mixed;

  }

  for ( int imix=0; imix<num_mixed; imix++ ) {

    double fraction = sdom->mixed_fraction[imix];

    if ( fraction == 1.0 ) {
	
      int zone = sdom->mixed_to_zones[imix];
      int material = sdom->mixed_material[imix];

      h_Aarray_scatter[num_zones_scatter] = d_inflated_material+num_moments*num_moments*material;
      h_Barray_scatter[num_zones_scatter] = sdom->d_phi+num_moments*num_groups*zone;
      h_Carray_scatter[num_zones_scatter] = sdom->d_phi_out+num_moments*num_groups*zone;
      
      num_zones_scatter++;

    }
    
  }

  cudaMemcpy(d_Aarray_scatter, h_Aarray_scatter, 3*num_mixed*sizeof(double*), cudaMemcpyHostToDevice);

  const double beta = 1.0;
  const double alpha = 1.0;
    
  cublasStatus_t custat;
  custat = cublasDgemmBatched ( kripke_cbhandle, CUBLAS_OP_N, CUBLAS_OP_N, num_moments, num_groups, num_moments, &alpha, 
				(const double **) d_Aarray_scatter, num_moments, 
				(const double **) d_Barray_scatter, num_moments, &beta, 
				d_Carray_scatter, num_moments, num_zones_scatter );
  
  if (custat) abort();
  cudaCheckError();


  // Then do all of the zones for which there is more than one material one at a time (slow).
  // (Common fractions could be pulled out, and batched - but seems it is less likely that a general model
  // would have common fractions?  A more permanent solution would be to just use a batched routine which
  // accepts a pointer to a list of alpha values - shouldn't be too hard.)
  for ( int imix=0; imix<num_mixed; imix++ ) {

    const double beta = 1.0;  

    double fraction = sdom->mixed_fraction[imix];
    
    if ( fraction != 1.0 ) {
      int zone = sdom->mixed_to_zones[imix];
      int material = sdom->mixed_material[imix];
      
      cublasStatus_t custat;
      
      custat = cublasDgemm_v2 ( kripke_cbhandle, CUBLAS_OP_N, CUBLAS_OP_N, num_moments, num_groups, num_moments, &fraction, 
				d_inflated_material+num_moments*num_moments*material, num_moments, 
				sdom->d_phi+num_moments*num_groups*zone, num_moments, &beta, 
				sdom->d_phi_out+num_moments*num_groups*zone, num_moments );
      if (custat) abort();
    }
  }
}



int cuda_LTimes_ZGD(double *d_phi, double *h_psi, double *d_psi, double *d_ell,
                    int num_zones, int num_groups, int num_local_directions, 
		    int num_local_groups, int nidx, int group0){

/*
 * This is just L * 'a subdomain worth of' Psi. 
 *
 *
 *
 *
 */

  cudaCheckError();
  int dim_y = 4;  
  dim3 threadsPerBlock(32,dim_y);

  cublasStatus_t custat;

  
  if ( ! d_psi ) {
    // d_psi gets created in sweep - so the first time this is called d_psi doesn't exist - fine just do it the old way

    LTimes_ZGD<<<num_zones,threadsPerBlock,nidx*dim_y*sizeof(double)>>>
      (d_phi,h_psi,d_ell,num_zones,num_groups,num_local_directions,num_local_groups,nidx,group0);
    
    cudaCheckError();
    
    custat = cublasCreate (&kripke_cbhandle);
    if ( custat ){
      printf ("Error cublasCreate \n");
      abort();
    }

  }

  else {

    // Loop over zones
    const double alpha = 1.0;
    const double beta =  1.0;

    static double ** h_Aarray = NULL;
    static double ** h_Barray = NULL;
    static double ** h_Carray = NULL;
    static double ** d_Aarray = NULL;
    static double ** d_Barray = NULL;
    static double ** d_Carray = NULL;
    static int ref_num_zones = 0;

    static cudaStream_t ltimesStream[8];
    static int ltimesstream = 0;

    if ( ! h_Aarray ) {

      // batched operation is the same every time - so only need to do this once.
      // small list of pointers to store on device
      // but means that full perf is not achieved until the 3rd iteration

      ref_num_zones = num_zones;
      h_Aarray = (double **) malloc (3*num_zones*sizeof(double*));
      h_Barray = h_Aarray + num_zones;
      h_Carray = h_Barray + num_zones; 
      
      cudaMalloc ((void**)&d_Aarray, 3*num_zones*sizeof(double*) ); 
      d_Barray = d_Aarray + num_zones;
      d_Carray = d_Barray + num_zones;
    }

    if ( num_zones != ref_num_zones ) {
      printf ("incorrect number of zones \n");
      abort();
    }

    // fix 12/26/16 - update the input lists every time
    // relative diff in results is now 1e-12.  as good as ever.

    for ( int izone=0; izone<num_zones; izone++ ) {
      h_Aarray[izone] = d_ell;
      h_Barray[izone] = d_psi+izone*num_local_groups*num_local_directions;
      h_Carray[izone] = d_phi+izone*num_groups*nidx+group0*nidx;
    }
      
    cudaMemcpy(d_Aarray, h_Aarray, 3*num_zones*sizeof(double*), cudaMemcpyHostToDevice);

    cudaError_t cuerr;
    for ( int i=0; i<8; i++ ) {
      cuerr = cudaStreamCreate ( &(ltimesStream[i]) );
      if ( cuerr ) abort();
    }

    cublasSetStream(kripke_cbhandle, ltimesStream[ltimesstream]);
    ltimesstream++;
    ltimesstream = ltimesstream%8;

    if ( 1 ) {
      custat = cublasDgemmBatched ( kripke_cbhandle, CUBLAS_OP_N, CUBLAS_OP_N, nidx, num_local_groups, 
				    num_local_directions, &alpha, (const double **) d_Aarray, nidx, 
				    (const double **) d_Barray, num_local_directions, &beta, 
				    d_Carray, nidx, num_zones );
      if ( custat ) {
	printf ("custat = %d \n", custat);
      }
    
    }
    else {

      static int * d_mlist = NULL;
      static int * d_nlist = NULL;
      static int * d_klist = NULL;
      static int * d_ldalist = NULL;
      static int * d_ldblist = NULL;
      static int * d_ldclist = NULL;

      if ( d_mlist == NULL ) {
	
	int *mlist, *nlist, *klist, *ldalist, *ldblist, *ldclist;
	mlist = (int *) malloc ( num_zones * sizeof(int));
	nlist = (int *) malloc ( num_zones * sizeof(int));
	klist = (int *) malloc ( num_zones * sizeof(int));
	ldalist = (int *) malloc ( num_zones * sizeof(int));
	ldblist = (int *) malloc ( num_zones * sizeof(int));
	ldclist = (int *) malloc ( num_zones * sizeof(int));
	
	cudaMalloc ((void**) &d_mlist, num_zones*sizeof(int));
	cudaMalloc ((void**) &d_nlist, num_zones*sizeof(int));
	cudaMalloc ((void**) &d_klist, num_zones*sizeof(int));
	cudaMalloc ((void**) &d_ldalist, num_zones*sizeof(int));
	cudaMalloc ((void**) &d_ldblist, num_zones*sizeof(int));
	cudaMalloc ((void**) &d_ldclist, num_zones*sizeof(int));
	
	printf (" m = %d, n = %d, k = %d \n", nidx, num_local_groups, num_local_directions);
	
	for ( int i=0; i<num_zones; i++ ) {
	  mlist[i] = nidx;
	  nlist[i] = num_local_groups;
	  klist[i] = num_local_directions;
	  ldalist[i] = nidx;
	  ldblist[i] = num_local_directions;
	  ldclist[i] = nidx;
	}
	
	cudaError_t cuerr8 = cudaGetLastError();
	if ( cuerr8 ) {
	  printf ("dangling error \n");
	}

	cudaMemcpy ( d_mlist, mlist, num_zones*sizeof(int), cudaMemcpyHostToDevice );
	cudaMemcpy ( d_nlist, nlist, num_zones*sizeof(int), cudaMemcpyHostToDevice );
	cudaMemcpy ( d_klist, klist, num_zones*sizeof(int), cudaMemcpyHostToDevice );
	cudaMemcpy ( d_ldalist, ldalist, num_zones*sizeof(int), cudaMemcpyHostToDevice );
	cudaMemcpy ( d_ldblist, ldblist, num_zones*sizeof(int), cudaMemcpyHostToDevice );
	cudaMemcpy ( d_ldclist, ldclist, num_zones*sizeof(int), cudaMemcpyHostToDevice );
	cudaDeviceSynchronize();
	cuerr8 = cudaGetLastError();
	if ( cuerr8 ) {
	  printf ("Error after cudaMemcpy \n");
	  abort();
	}

	free ( mlist );
	free ( nlist );
	free ( klist );
	free ( ldalist );
	free ( ldblist );
	free ( ldclist );

      }
      
      //printf ("Calling simple 1block batch \n");

      dgemm_custom_simple_1block_batch ( ltimesStream[ltimesstream],
					 CUBLAS_OP_N,
					 CUBLAS_OP_N,
					 d_mlist,
					 d_nlist,
					 d_klist,
					 &alpha,
					 (const double **) d_Aarray,
					 d_ldalist,
					 (const double **) d_Barray, 
					 d_ldblist,
					 &beta,
					 d_Carray, 
					 d_ldclist,
					 num_zones );

      //cudaDeviceSynchronize();
      cudaError_t cuerr8 = cudaGetLastError();
      if ( cuerr8 ) {
	printf ("Error after 1block batch \n");
	//abort();
      }
      
    }
    
    // terminate cuBLAS - need to move to some cleanup routine
    //cublasDestroy (cbhandle);
    
    // runtime += omp_get_wtime();
    // printf ("cublasDgemmBatche LTimes performance = %e GF/s\n", 2.0e-9*nidx*num_local_groups*num_local_directions*num_zones/runtime);
      
    cudaCheckError();

  }

  return 0;
}

__global__ void  LTimes_ZGD(double *phi, double * __restrict__ psi, const double * __restrict__ ell,
                                int num_zones, int num_groups, int num_local_directions, int num_local_groups, 
                                int nidx, int group0){


      extern __shared__ double ss_phi[];

      int z = blockIdx.x;
      double * __restrict__ block_phi = &phi[z*num_groups*nidx + group0*nidx];
      double * __restrict__ block_psi = &psi[z*num_local_groups*num_local_directions];

#if 1
     for(int group = threadIdx.y ;group < num_local_groups; group += blockDim.y){

        double *ss_phi_group = &ss_phi[nidx*threadIdx.y];


//   offset = z * zones*groups +
//            g * directions +
//            d;

        for(int nm_offset = threadIdx.x; nm_offset < nidx; nm_offset+=blockDim.x)
           ss_phi_group[nm_offset] = block_phi[nm_offset+nidx*group];

        for (int d = 0; d < num_local_directions; d++) {

          double psi_d = block_psi[d+group*num_local_directions];

          for(int nm_offset = threadIdx.x; nm_offset < nidx; nm_offset+=blockDim.x)
            ss_phi_group[nm_offset] += ell[nm_offset + nidx*d] * psi_d;
        }

        for(int nm_offset = threadIdx.x; nm_offset < nidx; nm_offset+=blockDim.x)
           block_phi[nm_offset+nidx*group] =  ss_phi_group[nm_offset];
      }

#else

      for(int group = 0;group < num_local_groups;++ group){

        for(int nm_offset = threadIdx.x; nm_offset < nidx; nm_offset+=blockDim.x)
           ss_phi[nm_offset] = block_phi[nm_offset+nidx*group];
 
        for (int d = 0; d < num_local_directions; d++) {

          double psi_d = block_psi[d+group*num_local_directions];
           
          for(int nm_offset = threadIdx.x; nm_offset < nidx; nm_offset+=blockDim.x)
            ss_phi[nm_offset] += ell[nm_offset + nidx*d] * psi_d;          
        }

        for(int nm_offset = threadIdx.x; nm_offset < nidx; nm_offset+=blockDim.x)
           block_phi[nm_offset+nidx*group] =  ss_phi[nm_offset];

      }

#endif

}


/*******************/
int  cuda_scattering_ZGD(int *d_mixed_to_zones, int *d_mixed_material, double *d_mixed_fraction, int *d_mixed_offset,
                      double *d_phi, double *d_phi_out, double *d_sigs0, double *d_sigs1, double *d_sigs2,
                      int *d_moment_to_coeff, int num_mixed, int num_moments, int num_groups, int num_coeff){


    int y_dim = 6;
    dim3 threadsPerBlock(32,y_dim);
    
    scattering_ZGD_step3<<<480,threadsPerBlock,num_groups*y_dim*sizeof(double)>>>(d_mixed_to_zones,d_mixed_material,d_mixed_fraction,d_mixed_offset,
										  d_phi,d_phi_out,d_sigs0, d_sigs1, d_sigs2, d_moment_to_coeff,num_mixed,num_moments,num_groups,num_coeff);

    cudaCheckError();
    
    return 0;
}

/*******************/



__global__ void scattering_ZGD(int * __restrict__ mixed_to_zones, int * __restrict__ mixed_material,
                               double * __restrict__ mixed_fraction, int * __restrict__ mixed_offset,
                               double * __restrict__ phi, double *phi_out, double * __restrict__ sigs0,
                               double * __restrict__ sigs1, double * __restrict__ sigs2,
                               int * __restrict__ moment_to_coeff, int num_mixed, int num_moments, int num_groups, int num_coeff){


   int mix_min = mixed_offset[blockIdx.x];
   int mix_max = mixed_offset[blockIdx.x+1];
   double *d_sigs[3];
   d_sigs[0] = sigs0;
   d_sigs[1] = sigs1;
   d_sigs[2] = sigs2;

   for(int mix = mix_min;mix < mix_max; ++mix){
      int zone = mixed_to_zones[mix];
      int material = mixed_material[mix];
      double fraction = mixed_fraction[mix];
      double *sigs_g_gp = d_sigs[material];
      double *phi_z_g = &phi[zone*num_groups*num_moments];
      double *phi_out_z_gp = &phi_out[zone*num_groups*num_moments];

//LG unroll the outer loop
      for(int g = 0; g < num_groups;++g){

        for(int gp = threadIdx.y; gp < num_groups; gp += blockDim.y){

          for(int nm = threadIdx.x; nm < num_moments; nm += blockDim.x){
            // map nm to n
            int n = moment_to_coeff[nm];

            phi_out_z_gp[nm+gp*num_moments] += 
                sigs_g_gp[n + gp*num_coeff + g*num_groups*num_coeff] * 
                phi_z_g[nm + g*num_moments] * fraction;
          }
        }
      }

   }
}



__global__ void scattering_ZGD_step2(const int * __restrict__ mixed_to_zones, const int * __restrict__ mixed_material,
                               const double * __restrict__ mixed_fraction, const int * __restrict__ mixed_offset,
                               const double * __restrict__ phi, double *phi_out, double * __restrict__ sigs0,
                               const double * __restrict__ sigs1, const double * __restrict__ sigs2,
                               const int * __restrict__ moment_to_coeff, int num_mixed, int num_moments, int num_groups, int num_coeff){

   extern __shared__ double phi_z_g_ss[];

   int mix_min = mixed_offset[blockIdx.x];
   int mix_max = mixed_offset[blockIdx.x+1];
   int gid = threadIdx.x + threadIdx.y*blockDim.x;
   const double *d_sigs[3];
   d_sigs[0] = sigs0;
   d_sigs[1] = sigs1;
   d_sigs[2] = sigs2;

   for(int mix = mix_min;mix < mix_max; ++mix){
      int zone = mixed_to_zones[mix];
      int material = mixed_material[mix];
      double fraction = mixed_fraction[mix];
      const double * __restrict__ sigs_g_gp = d_sigs[material];
      const double * __restrict__ phi_z_g = &phi[zone*num_groups*num_moments];
      double * __restrict__ phi_out_z_gp = &phi_out[zone*num_groups*num_moments];

//LG unroll the outer loop
      for(int g = 0; g < num_groups;++g){

        __syncthreads();
        for(int nm = gid; nm < num_moments; nm += blockDim.x*blockDim.y)
          phi_z_g_ss[nm] =  phi_z_g[nm + g*num_moments]*fraction;
         __syncthreads(); 

        for(int gp = threadIdx.y; gp < num_groups; gp += blockDim.y){

          for(int nm = threadIdx.x; nm < num_moments; nm += blockDim.x){
            // map nm to n
            const int n = moment_to_coeff[nm];

            phi_out_z_gp[nm+gp*num_moments] +=
                sigs_g_gp[n + gp*num_coeff + g*num_groups*num_coeff] *
                phi_z_g_ss[nm];
          }
        }
      }

   }
}


__global__ void scattering_ZGD_step3(const int * __restrict__ mixed_to_zones, const int * __restrict__ mixed_material,
                               const double * __restrict__ mixed_fraction, const int * __restrict__ mixed_offset,
                               const double * __restrict__ phi, double *phi_out, double * __restrict__ sigs0,
                               const double * __restrict__ sigs1, const double * __restrict__ sigs2,
                               const int * __restrict__ moment_to_coeff, int num_mixed, int num_moments, int num_groups, int num_coeff){

   extern __shared__ double phi_out_ss[];

   int mix_min = mixed_offset[blockIdx.x];
   int mix_max = mixed_offset[blockIdx.x+1];
   const double *d_sigs[3];
   d_sigs[0] = sigs0;
   d_sigs[1] = sigs1;
   d_sigs[2] = sigs2;

   for(int mix = mix_min;mix < mix_max; ++mix){

      int zone = mixed_to_zones[mix];
      int material = mixed_material[mix];
      double fraction = mixed_fraction[mix];
      const double * __restrict__ sigs_g_gp = d_sigs[material];
      const double * __restrict__ phi_z_g = &phi[zone*num_groups*num_moments];
      double * __restrict__ phi_out_z_gp = &phi_out[zone*num_groups*num_moments];

//LG unroll the outer loop

      for(int gp = threadIdx.y; gp < num_groups; gp += blockDim.y){

        double *phi_out_ss_gp = &phi_out_ss[num_groups*threadIdx.y];
        for(int nm = threadIdx.x; nm < num_moments; nm += blockDim.x)
           phi_out_ss_gp[nm] =  phi_out_z_gp[nm+gp*num_moments];

        for(int g = 0; g < num_groups;++g){

          const int nm_shift = g*num_groups*num_coeff + gp*num_coeff ;

          for(int nm = threadIdx.x; nm < num_moments; nm += blockDim.x){
            // map nm to n
            const int n = moment_to_coeff[nm];

            phi_out_ss_gp[nm] +=
                sigs_g_gp[n + nm_shift] *
                phi_z_g[nm + g*num_moments]*fraction;

          }
        }
        for(int nm = threadIdx.x; nm < num_moments; nm += blockDim.x)
           phi_out_z_gp[nm+gp*num_moments] = phi_out_ss_gp[nm];

      }

   }
}







__global__ void source_ZGD2 ( int * __restrict__ mixed_to_zones, int * __restrict__ mixed_material,
			      double * __restrict__ mixed_fraction, /* int * __restrict__ mixed_offset, */
			      double *phi_out, int num_moments, int num_groups, const int num_mixed )
/*
  Seems simpler to understand thatn source_ZGD.
  Doesn't used mixed_offset - which was hard to follow.  
  One thread per zone.
  Each thread loops over num_groups - and updates the zero moment for each group.
  Inner dimension is moments - so no possibility of coalesced memory access.
  But this is very fast - ~200 usecs.  Not important at all (at the time of this writing).
*/
{


  int tid = threadIdx.x + blockIdx.x*blockDim.x;

  if ( tid < num_mixed ) {
  
    int material = mixed_material[tid];

    if(material == 0){

      int zone = mixed_to_zones[tid];
      double fraction = mixed_fraction[tid];
      double *phi_out_z = &phi_out[zone*num_moments*num_groups];

      for ( int g=0; g<num_groups; g++ ){
	phi_out_z[g*num_moments] += 1.0 * fraction;
      }

    }
  }
}


int  cuda_source_ZGD(int *d_mixed_to_zones, int *d_mixed_material, double *d_mixed_fraction, int *d_mixed_offset,
		     double *d_phi_out, int num_moments, int num_groups, int num_mixed )
{


  if ( 0 ) {
    
    dim3 threadsPerBlock(32,1);
    source_ZGD<<<480,threadsPerBlock>>>(d_mixed_to_zones,d_mixed_material,d_mixed_fraction,d_mixed_offset,
					d_phi_out,num_moments,num_groups);
    
    cudaCheckError();
  }
  else {
    int nthreads = 256;
    dim3 blocks( (num_mixed+nthreads-1)/nthreads, 1, 1);
    source_ZGD2 <<<blocks,nthreads>>> ( d_mixed_to_zones, d_mixed_material, d_mixed_fraction, 
					d_phi_out, num_moments, num_groups, num_mixed);
  }
  return 0;
}

__global__ void source_ZGD(int * __restrict__ mixed_to_zones, int * __restrict__ mixed_material,
                               double * __restrict__ mixed_fraction, int * __restrict__ mixed_offset,
                               double *phi_out, int num_moments, int num_groups){


   int mix_min = mixed_offset[blockIdx.x];
   int mix_max = mixed_offset[blockIdx.x+1];


   for(int mix = mix_min;mix < mix_max; ++mix){
  
      int material = mixed_material[mix];

      if(material == 0){
          int zone = mixed_to_zones[mix];
          double fraction = mixed_fraction[mix];
          double *phi_out_z = &phi_out[zone*num_moments*num_groups];
          for(int g = threadIdx.x; g < num_groups; g += blockDim.x){
            phi_out_z[g*num_moments] += 1.0 * fraction;
          }
       }
    }

}


int  cuda_LPlusTimes_ZGD(double *d_rhs, double *d_phi_out, double *d_ell_plus,
                    int num_zones, int num_groups, int num_local_directions, int num_local_groups,
                    int nidx, int group0){


#ifdef CU_TIMING
  cudaEvent_t start,stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  float time_ms, time_s;
  cudaEventRecord(start);
#endif
  
  cublasStatus_t custat;


#ifndef KRIPKE_ZGD_FLUX_REGISTERS
    
    dim3 threadsPerBlock(32);

    LPlusTimes_ZGD<<<num_zones,threadsPerBlock,num_local_directions*sizeof(double)>>>(d_rhs,d_phi_out,d_ell_plus,num_zones,num_groups, 
										      num_local_directions, num_local_groups, nidx, group0);

#else
    
    const double alpha = 1.0;
    const double beta =  1.0;

    static double ** h_Aarray_lpt = NULL;
    static double ** h_Barray_lpt = NULL;
    static double ** h_Carray_lpt = NULL;
    static double ** d_Aarray_lpt = NULL;
    static double ** d_Barray_lpt = NULL;
    static double ** d_Carray_lpt = NULL;
    static int ref_num_zones_lpt = 0;

    if ( ! h_Aarray_lpt ) {

      // batched operation is the same every time - so only need to do this once.
      // small list of pointers to store on device
      // but means that full perf is not achieved until the 3rd iteration

      ref_num_zones_lpt = num_zones;
      h_Aarray_lpt = (double **) malloc (3*num_zones*sizeof(double*));
      h_Barray_lpt = h_Aarray_lpt + num_zones;
      h_Carray_lpt = h_Barray_lpt + num_zones;
      
      cudaMalloc ((void**)&d_Aarray_lpt, 3*num_zones*sizeof(double*) ); 
      d_Barray_lpt = d_Aarray_lpt + num_zones;
      d_Carray_lpt = d_Barray_lpt + num_zones;

    }
    
    if ( num_zones != ref_num_zones_lpt ) {
      printf ("incorrect number of zones in LPlusTimes \n");
      abort();
    }
    
    for ( int izone=0; izone<num_zones; izone++ ) {
      h_Aarray_lpt[izone] = d_ell_plus;
      h_Barray_lpt[izone] = d_phi_out + group0*nidx + izone*num_groups*nidx;
      h_Carray_lpt[izone] = d_rhs+izone*num_local_groups*num_local_directions;
    }
      
    cudaMemcpy(d_Aarray_lpt, h_Aarray_lpt, 3*num_zones*sizeof(double*), cudaMemcpyHostToDevice);

    custat = cublasDgemmBatched ( kripke_cbhandle, CUBLAS_OP_T, CUBLAS_OP_N, num_local_directions, num_local_groups, nidx, &alpha, (const double **) d_Aarray_lpt, nidx, 
    				  (const double **) d_Barray_lpt, nidx, &beta, 
    				  d_Carray_lpt, num_local_directions, num_zones );
    
    if ( custat ) {
      printf ("custat = %d \n", custat);
    }
    
    // terminate cuBLAS - need to move to some cleanup routine
    //cublasDestroy (cbhandle);
    
    cudaDeviceSynchronize();

    cudaCheckError();

#endif

  #ifdef CU_TIMING
  cudaEventRecord(stop);
  cudaDeviceSynchronize();
  cudaCheckError();
  cudaEventElapsedTime(&time_ms,start,stop);
  time_s=time_ms*.001;
  printf("ZGD: time to LPlusTimes_ZGD (GPU): %g [s]\n",time_s);
  #endif

  return 0;

}



__global__ void  LPlusTimes_ZGD(double *rhs, double * __restrict__ phi_out,
                                double * __restrict__ ell_plus,
                                int num_zones, int num_groups, int num_local_directions, int num_local_groups, 
                                int nidx, int group0){


//LG consider running  for(int nm_offset  in parallel and using reduction within a warp

      int z = blockIdx.x;
      extern __shared__ double rhs_acc[];
      double *block_rhs =  &rhs[z*num_local_groups*num_local_directions];


//   offset = z * zones*groups +
//            g * directions +
//            d;


      double *block_phi_out = &phi_out[z*num_groups*nidx + group0*nidx];

      for(int group = 0; group < num_local_groups;++ group){

        for (int d = threadIdx.x; d < num_local_directions; d+=blockDim.x) {
          rhs_acc[d] = 0.0;

          for(int nm_offset = 0;nm_offset < nidx;++nm_offset)
             rhs_acc[d] += ell_plus[nm_offset + d*nidx] * block_phi_out[nm_offset + group*nidx];
          
          block_rhs[d+num_local_directions*group] += rhs_acc[d];
        }
      }
}




template <int iOct, int jOct, int kOct>
__global__ void sweep_over_hyperplane_ZGD_fluxRegisters ( const int nBlocks_j,
							  const int nBlocks_k,
							  const int i_inc,
							  const int j_inc,
							  const int k_inc,
							  const int direction_offsetX,
							  const int groupX, 
							  const int num_groups,
							  const int num_directions,
							  const int local_imax,
							  const int local_jmax,
							  const int local_kmax,
							  double * __restrict__ d_dx, 
							  double * __restrict__ d_dy, 
							  double * __restrict__ d_dz, 
							  double * __restrict__ d_rhs, 
							  const double * __restrict__ d_sigt, 
							  Directions * __restrict__  d_direction,
							  double * __restrict__ d_psi, 
							  double * __restrict__ flux_boundary_i,
							  double * __restrict__ flux_boundary_j,
							  double * __restrict__ flux_boundary_k,
							  int kernelIn, 
							  int nKernels
							  )

/*

  Each block will process 16 directions per zone for FluxPlaneDim x FluxPlaneDim x local_imax zones 

  DRAM traffic is minimized by keeping x fluxes local to the thread (registers) and 
  sharing y and z fluxes using SMEM.

  smem_flux_j = 8 x 8 x 16 * 8 bytes = 8k Bytes arranged as dir / j / k  (say, for 8, 16)
  smem_flux_k = 8 x 8 x 16 * 8 bytes = 8k Bytes arranged as dir / j / k 

  Issues:
      1. poor hit-rate for sigma - due to it being ordered fastest by groups (and a kernel only does one group)
      2. need to advance-load flux_boundary_j and _k.  Some loads per hyperplane slow down almost all hyperplanes.
      3. ? maybe advance-load sigma?
      4. double-buffer j/k-plane data.  This would remove the need for a syncthreads (but occupy double the SMEM)

*/

{

  int jBlock, kBlock;
  int i, z;
  double flux_i, flux_j, flux_k;

  //	      (kernel%(num_directions/FluxGroupDim))*FluxGroupDim,
  //	      kernel/(num_directions/FluxGroupDim),

  
  int kernel, direction_offset, group;

  kernel = kernelIn + blockIdx.x;

  if ( kernel >= nKernels ) {
    return;
  }

  direction_offset = (kernel%(num_directions/FluxGroupDim))*FluxGroupDim;
  group = kernel/(num_directions/FluxGroupDim);
  

  // Octant-specific rules
  if ( iOct == 0 ) {
    i = 0;
  }
  else if ( iOct == 1 ) {
    i = local_imax-1;
  }

  if ( jOct == 0 ) {
    jBlock = 0;
  }
  else {
    jBlock = nBlocks_j - 1;
  }

  if ( kOct == 0 ) {
    kBlock = 0;
  }
  else {
    kBlock = nBlocks_k - 1;
  }

  // setup space to exchange j, k fluxes in SMEM  - per block
  extern __shared__ double smem[];
 
  double * smem_flux_j = (double*) smem;
  double * smem_flux_k = (double*) &smem_flux_j [FluxPlaneDim*FluxPlaneDim*FluxGroupDim];


  const int tid = threadIdx.x ;                                               // local (i.e. rank) thread index
  const int d = tid%FluxGroupDim + direction_offset;                          // direction index (1/2 warp)

  // block index (i.e. within local_imax x 8 x 8 block)
  const int dd = tid%FluxGroupDim;                                                       // direction index (1/2 warp)
  const int jj = (tid/FluxGroupDim)%FluxPlaneDim;                                        // j location within 8 x 8 plane
  const int kk = (tid/FluxGroupDim)/FluxPlaneDim;                                        // k location within 8 x 8 plane

  // node index (i.e. within local_imax x local_jmax x local_kmax)
  int j = jBlock*FluxPlaneDim+jj;                                              // local (i.e. rank) y zone
  int k = kBlock*FluxPlaneDim+kk;                                              // local (i.e. rank) z zone

  // load this thread's constant data                                          // i.e. constant for this zone-pencil in x
  const double xcos_dxi = d_direction [d].xcos * 2.0 / d_dx[i+1];              // zero effect on performance
  const double ycos_dyj = d_direction [d].ycos * 2.0 / d_dy[j+1];
  const double zcos_dzk = d_direction [d].zcos * 2.0 / d_dz[k+1];
  const int dir_grp = num_directions * num_groups;
  const int gd = d + group * num_directions;

  const double * __restrict__ block_sigt = NULL;  
  double * __restrict__ block_rhs = NULL;    // pointer to rhs data
  double * __restrict__ block_psi;           // pointer to psi data

  // load in the incoming i-plane flux from GMEM for each thread
  flux_i = flux_boundary_i [dir_grp*(j+k*local_jmax) + group*num_directions + d ];

  // initialize the pointers to GMEM data
  z = j*local_imax + k*local_imax*local_jmax + i;
  block_sigt = &d_sigt[z*num_groups+group];
  block_rhs = &d_rhs[z*dir_grp + gd];    // pointer to rhs data
  block_psi = &d_psi[z*dir_grp + gd];    // pointer to psi data

  // loop over i-planes 
  // (a.k.a. loop over all the hyperplanes in the block allocated to the NODE)
  // This is done to minimize 'tail' effects. 8x8 pencils stack up with each other and leave zero gaps. 
  // Applies to all 8x8 pencils in the domain for this set of 16 directions.
  // A 32 x 32 node domain has 16 8x8 pencils.  Scaned linearly by j then by k.
  // Certainly, finding ways to parallelize this across multiple blocks would be a further optimization.
  //for ( int hplane=0; hplane < local_imax*nBlocks_j*nBlocks_k + FluxGroupDim; ++hplane ) {
  for ( int hplane=0; hplane < local_imax*nBlocks_j*nBlocks_k + 2*FluxPlaneDim; ++hplane ) {

    // master sync - basically between hyperplanes - required to ensure that all flux data is in SMEM.
    // Significant performance limiter.  Removing this gives a 20% performance boost.
    __syncthreads();

    int hptest = 0;
    if ( jOct == 0 ) {
      hptest += jj;
    }
    else {
      hptest += (FluxPlaneDim-1-jj);
    }

    if ( kOct == 0 ) {
      hptest += kk;
    }
    else {
      hptest += (FluxPlaneDim-1-kk);
    }

    // check to see if the current zone is in the current hyperplane
    if ( kBlock >= 0 && kBlock < nBlocks_k && hplane >= hptest ) {      // this only affects the start and stop.  they start at hplane >= jj+kk
                                                                        // once the threads have been started, they are active until they stop at kBlock == nBlocks_k
                                                                        // synchthreads keeps the started hyperplane in step

      if ( jOct == 0 ) {
	if ( jj == 0 ) {
	  // on the j-plane input boundary, so get flux_j from GMEM
	  flux_j = flux_boundary_j [dir_grp*(i+k*local_imax) + group*num_directions + d]; 
	} 
	else {
	  // get flux_j from one block to the side in SMEM
	  flux_j = smem_flux_j[FluxGroupDim*(kk+jj*FluxPlaneDim - FluxPlaneDim)+dd];
	}
      }
      else {
	if ( jj == FluxPlaneDim-1 ) {
	  // on the j-plane input boundary, so get flux_j from GMEM
	  flux_j = flux_boundary_j [dir_grp*(i+k*local_imax) + group*num_directions + d]; 
	} 
	else {
	  // get flux_j from one block to the side in SMEM
	  flux_j = smem_flux_j[FluxGroupDim*(kk+jj*FluxPlaneDim + FluxPlaneDim)+dd];
	}
      }

      if ( kOct == 0 ) {
	if ( kk == 0 ) {
	  // on the k-plane input boundary, so get flux_k from GMEM
	  flux_k = flux_boundary_k [dir_grp*(i+j*local_imax) + group*num_directions + d]; 
	}
	else {
	  // get flux_k from the block in the lower row in SMEM
	  flux_k = smem_flux_k [FluxGroupDim*(jj+kk*FluxPlaneDim - FluxPlaneDim)+dd];
	}
      }
      else {
	if ( kk == FluxPlaneDim-1 ) {
	  // on the k-plane input boundary, so get flux_k from GMEM
	  flux_k = flux_boundary_k [dir_grp*(i+j*local_imax) + group*num_directions + d]; 
	}
	else {
	  // get flux_k from the block in the lower row in SMEM
	  flux_k = smem_flux_k [FluxGroupDim*(jj+kk*FluxPlaneDim + FluxPlaneDim)+dd];
	}
      }

      // calculate the new zonal flux
      double psi_z_g_d = ( ( 
			    *block_rhs
			    //__ldg(block_rhs)
			    + flux_i * xcos_dxi
			    + flux_j * ycos_dyj
			    + flux_k * zcos_dzk ) / 
			    ( 
			    *block_sigt  
			    //__ldg(block_sigt) 
			    + xcos_dxi 
			    + ycos_dyj  
			    + zcos_dzk )  
			   );

      // output new Psi to GMEM
      *block_psi = psi_z_g_d;

      psi_z_g_d *= 2;

      // update flux-i in register only
      flux_i = psi_z_g_d - flux_i;

      // update flux-j,k   
      flux_j = psi_z_g_d - flux_j;
      flux_k = psi_z_g_d - flux_k; 
    }

    __syncthreads();        // make sure all reads have completed before beginning writes

    if ( kBlock >= 0 && kBlock < nBlocks_k && hplane >= hptest ) {       

      // output flux boundaries
      if ( jOct == 0 ) {
	if ( jj == FluxPlaneDim-1 ) {
	  // on the j-plane output boundary, so send flux_j to GMEM 
	  flux_boundary_j [dir_grp*(i+k*local_imax) + group*num_directions + d] = flux_j; 
	}
	else {
	  // send flux_j to SMEM
	  if ( ( FluxGroupDim* ( kk+jj*FluxPlaneDim ) + dd < FluxPlaneDim*FluxPlaneDim*FluxGroupDim ) && (FluxGroupDim*(kk+jj*FluxPlaneDim)+dd  >= 0) ) {
	    smem_flux_j[FluxGroupDim*(kk+jj*FluxPlaneDim)+dd] = flux_j;	    
	  }	    
	}
      }
      else {
	if ( jj == 0 ) {
	  // on the j-plane output boundary, so send flux_j to GMEM
	  flux_boundary_j [dir_grp*(i+k*local_imax) + group*num_directions + d] = flux_j; 
	}
	else {
	  // send flux_j to SMEM
	  if ( (FluxGroupDim*(kk+jj*FluxPlaneDim)+dd) >= 0 &&  (FluxGroupDim*(kk+jj*FluxPlaneDim)+dd) < FluxPlaneDim*FluxPlaneDim*FluxGroupDim ) {
	    smem_flux_j[FluxGroupDim*(kk+jj*FluxPlaneDim)+dd] = flux_j;
	  }
	}
      }

      if ( kOct == 0 ) {
	if ( kk == FluxPlaneDim-1 ) {
	  // on the k-plane output boundary, so send flux_k to GMEM
	  flux_boundary_k [dir_grp*(i+j*local_imax) + group*num_directions + d] = flux_k; 
	}
	else {
	  // send flux_k to SMEM
	  if ( (FluxGroupDim*(jj+kk*FluxPlaneDim)+dd) >= 0 &&  (FluxGroupDim*(jj+kk*FluxPlaneDim)+dd) < FluxPlaneDim*FluxPlaneDim*FluxGroupDim ) {
	    smem_flux_k[FluxGroupDim*(jj+kk*FluxPlaneDim)+dd] = flux_k;
	  }
	}
      }
      else {
	if ( kk == 0 ) {
	  // on the k-plane output boundary, so send flux_k to GMEM
	  flux_boundary_k [dir_grp*(i+j*local_imax) + group*num_directions + d] = flux_k; 
	}
	else {
	  // send flux_k to SMEM	  
	  if ( (FluxGroupDim*(jj+kk*FluxPlaneDim)+dd) >= 0 && (FluxGroupDim*(jj+kk*FluxPlaneDim)+dd) < FluxPlaneDim*FluxPlaneDim*FluxGroupDim ) {
	    smem_flux_k[FluxGroupDim*(jj+kk*FluxPlaneDim)+dd] = flux_k;
	  }
	}
      }

      if ( iOct == 0 ) {//---------------------------------------------------------------------------------------------------------
	// translate down the x-direction by 1
	i++;
	
	// update zone for input (rhs) and output (psi)
	block_rhs += dir_grp;
	block_psi += dir_grp;   
	block_sigt += num_groups;

	// handle reaching the end of the i-domain
	if ( i == local_imax ) {
	  
	  // at end of i-domain, send flux_i to GMEM
	  flux_boundary_i [dir_grp*(j+k*local_jmax) + d + num_directions*group] = flux_i;  
	  
	  i = 0;                                           // i is reset

 	  if ( jOct == 0 ) {
	    j += FluxPlaneDim;                                          // increment to the same jj in the next j-block
	    jBlock++;                                       
	    
	    // handle reaching the end of the j-domain
	    if ( j >= local_jmax ) {
	      j = jj;
	      jBlock = 0;
	      if ( kOct == 0 ) {
		k += FluxPlaneDim;
		kBlock++;
	      }
	      else {
		k -= FluxPlaneDim;
		kBlock--;
	      }
	    }
	  }
	  else {
	    j -= FluxPlaneDim;                                          // increment to the same jj in the next j-block
	    jBlock--;                                       
	    
	    // handle reaching the end of the j-domain
	    if ( j < 0 ) {
	      jBlock = nBlocks_j-1;
	      j = jBlock*FluxPlaneDim+jj;
	      if ( kOct == 0 ) {
		k += FluxPlaneDim;
		kBlock++;
	      }
	      else {
		k -= FluxPlaneDim;
		kBlock--;
	      }
	    }
	  }

	  // update rhs, psi, sigt pointers to GMEM
	  z = j*local_imax + k*local_imax*local_jmax + i;
	  block_sigt = &d_sigt[z*num_groups+group];
	  block_rhs = &d_rhs[z*dir_grp + gd];              // pointer to rhs data
	  block_psi = &d_psi[z*dir_grp + gd];              // pointer to psi data
	  
	  // load the input flux_i for the new 8x8 block.
	  flux_i = flux_boundary_i [dir_grp*(j+k*local_jmax) + group*num_directions + d ];
	
	}

      }

      else if ( iOct == 1 ) {//---------------------------------------------------------------------------------------------------------
	// translate down the x-direction by 1
	i--;
	
	// update zone for input (rhs) and output (psi)
	block_rhs -= dir_grp;
	block_psi -= dir_grp;   
	block_sigt -= num_groups;

	// handle reaching the end of the i-domain
	if ( i < 0 ) {
	  
	  // at end of i-domain, send flux_i to GMEM
	  flux_boundary_i [dir_grp*(j+k*local_jmax) + d + num_directions*group] = flux_i;  
	  
	  i = local_imax-1;                                // i is reset

	  if ( jOct == 0 ) {
	    j += FluxPlaneDim;                                          // increment to the same jj in the next j-block
	    jBlock++;                                       
	    
	    // handle reaching the end of the j-domain
	    if ( j >= local_jmax ) {
	      j = jj;
	      jBlock = 0;
	      if ( kOct == 0 ) {
		k += FluxPlaneDim;
		kBlock++;
	      }
	      else {
		k -= FluxPlaneDim;
		kBlock--;
	      }
	    }
	  }
	  else {
	    j -= FluxPlaneDim;                                          // increment to the same jj in the next j-block
	    jBlock--;                                       
	    
	    // handle reaching the end of the j-domain
	    if ( j < 0 ) {
	      jBlock = nBlocks_j-1;
	      j = jBlock*FluxPlaneDim+jj;
	      if ( kOct == 0 ) {
		k += FluxPlaneDim;
		kBlock++;
	      }
	      else {
		k -= FluxPlaneDim;
		kBlock--;
	      }
	    }
	  }

	  // update rhs, psi, sigt pointers to GMEM
	  z = j*local_imax + k*local_imax*local_jmax + i;
	  block_sigt = &d_sigt[z*num_groups+group];
	  block_rhs = &d_rhs[z*dir_grp + gd];              // pointer to rhs data
	  block_psi = &d_psi[z*dir_grp + gd];              // pointer to psi data
	  
	  // load the input flux_i for the new 8x8 block.
	  flux_i = flux_boundary_i [dir_grp*(j+k*local_jmax) + group*num_directions + d ];
	
	}

      }

    }

  }

}



int cuda_sweep_ZGD_fluxRegisters ( const int local_imax,
				   const int local_jmax,
				   const int local_kmax,
				   const int num_zones,
				   const int num_directions,
				   const int num_groups,
				   double * __restrict__ d_rhs,
				   const double * __restrict__ d_sigt,
				   Directions * __restrict__ d_direction,
				   double * __restrict__ d_dx,
				   double * __restrict__ d_dy,
				   double * __restrict__ d_dz,
				   double * __restrict__ h_psi,
				   double * __restrict__ h_i_plane,
				   double * __restrict__ h_j_plane,
				   double * __restrict__ h_k_plane,
				   int i_inc,
				   int j_inc,
				   int k_inc,
				   Subdomain * __restrict__ sdom
				   )

/* 

   Perform the sweep over the local zones while keeping the fluxes in register/smem.

   This removes the majority of the DRAM reads/writes (removes 6 of 8[9]) and should
   give a corresponding performance improvement due to sweep being entirely bandwidth bound.

   This sweep occurs on a single GPU.  Expect this to be managed by a single rank.

   *** Requirements: *** 
      #directions must be a multiple of FluxGroupDim
      #zones in y and z must be a multiple of FluxPlaneDim

   Each BLOCK handles FluxPlaneDim x FluxPlaneDim (j,k - zones) x FluxGroupDim (directions) x 1 (energy group) [need to generalize directions/groups]

   One CUDA Block loops serially through j-k 'pencils' for FluxGroupDim directions.
  
   Launch ngroups * ndirs / FluxGroupDim number of KERNELS.
   
   Keep x-flux in register.

   Communicate y and z fluxes, within a CUDA block, via smem.  
      16k Bytes per block 

   Synchronize within a block using syncthreads();

   Synchronize between MPI ranks as current.

*/

{

  cudaCheckError();
  // Input Checks - since the current implementation has some restrictions on input parameters
  {
#ifdef USE_PSI_HOST_MEM
    printf ("calling cuda_sweep_ZGD_fluxRegisters with USE_PSI_HOST_MEM not supported (sensible).\n");
    abort();
#endif

#ifdef USE_IJK_PLANE_HOST_MEM
    printf ("calling cuda_sweep_ZGD_fluxRegisters with USE_IJK_PLANE_HOST_MEM not supported (sensible).\n");
    abort();
#endif

    if ( local_jmax%FluxPlaneDim || local_kmax%FluxPlaneDim ) {
      printf ("local y and z zone extents must be multiples of 8: %d, %d \n", local_jmax, local_kmax);
      abort();
    }
    
    if ( num_directions%FluxGroupDim ) {
      printf ("number of directions MUST be a multiple of %d \n", FluxGroupDim);
      abort();
    }
  }

  // Allocate space on the GPU for Psi and copy Psi to GPU.
  double *d_psi = sdom->d_psi;
  size_t N = num_zones * num_directions * num_groups;

  {

#ifdef CU_TIMING
    cudaEventRecord(start);
#endif

#ifdef CU_TIMING
    cudaEventRecord(stop);
    cudaDeviceSynchronize();
    cudaCheckError();
    cudaEventElapsedTime(&time_ms,start,stop);
    time_s=time_ms*.001;
    printf("ZGD: time to copy PSI H2D: %g [s]\n",time_s);
#endif

  }

  // Allocate space on the GPU for input fluxes.
  // *** NOTE that this had previously been arranged as ZDG.  However, it is advantageous
  //          to re-arrange as ZGD since we will always be reading 16 directions at once and
  //          this ZGD arrangements results in completely coalesced loads.  BUT this adds a
  //          lot of overhead to the CPU-portion which would need to be removed to get 
  //          legitimate performance numbers at scale.

  size_t groups_dirs = num_directions * num_groups;
  int i_plane_zones = local_jmax * local_kmax * groups_dirs;
  int j_plane_zones = local_imax * local_kmax * groups_dirs;
  int k_plane_zones = local_imax * local_jmax * groups_dirs; 

  double *d_i_plane,  *d_j_plane, *d_k_plane;

  #ifdef CU_TIMING
  cudaEventRecord(start);
  #endif

  d_i_plane = sdom->d_i_plane;
  d_j_plane = sdom->d_j_plane;
  d_k_plane = sdom->d_k_plane;

  //double * tempi = sdom->tempi;
  //double * tempj = sdom->tempj;
  //double * tempk = sdom->tempk;

  cudaCheckError();

  //cudaEventRecord( sdom->sweepEvents[0], sdom->sweepStreams[0] );

  //for ( int i=0; i<32; i++ ) {
  //  cudaStreamWaitEvent( sdom->sweepStreams[i], sdom->sweepEvents[0], 0 );
  //}

  cudaCheckError();

  #ifdef CU_TIMING
  cudaEventRecord(stop);
  cudaDeviceSynchronize();
  cudaCheckError();
  cudaEventElapsedTime(&time_ms,start,stop);
  time_s=time_ms*.001;
  printf("ZGD: time to copy ijk_PLANE H2D: %g [s]\n",time_s);
  #endif

  // number of required blocks
  int nBlocks_j = local_jmax/FluxPlaneDim;
  int nBlocks_k = local_kmax/FluxPlaneDim;

  // number of required kernels
  size_t nKernels = num_groups * (num_directions / FluxGroupDim);

  int threadsPerBlock = FluxPlaneDim*FluxPlaneDim*FluxGroupDim;                                                     // !fixed! 8 x 8 x 16


  int octant = 0;
  if ( i_inc ) octant += 100;
  if ( j_inc ) octant += 10;
  if ( k_inc ) octant += 1;

  if ( 1 ) { 

  double mpi_time_start = MPI_Wtime();

  // Required shared memory to hold y and z fluxes (* group size * bytes).  
  int required_smem = FluxPlaneDim * FluxPlaneDim * FluxGroupDim * 8 * 2 ;

  double omptime = -1.0 * omp_get_wtime();

  cudaCheckError();

  for ( int kernel=0; kernel<nKernels; kernel+=NumSMs ) {

    if ( i_inc == 0 ) {
      if ( j_inc == 0 ) {
	if ( k_inc == 0 ) {
	  if ( kernel%NumSMs == 0 ) {
  cudaCheckError();
	    sweep_over_hyperplane_ZGD_fluxRegisters <0,0,0>
	      //	      <<< NumSMs, threadsPerBlock, required_smem, sdom->sweepStreams[kernel%32] >>> 
	      <<< NumSMs, threadsPerBlock, required_smem, sdom->subDStream >>> 
	      ( nBlocks_j, nBlocks_k, i_inc, j_inc, k_inc, (kernel%(num_directions/FluxGroupDim))*FluxGroupDim, 
		kernel/(num_directions/FluxGroupDim), num_groups, num_directions,
		local_imax, local_jmax, local_kmax, d_dx, d_dy, d_dz, d_rhs, d_sigt, d_direction, d_psi, d_i_plane, d_j_plane, d_k_plane, kernel, nKernels);
  cudaCheckError();
	  }
	}
	else {
	  if ( kernel%NumSMs == 0 ) {
  cudaCheckError();
	  sweep_over_hyperplane_ZGD_fluxRegisters <0,0,1>
	    <<< NumSMs, threadsPerBlock, required_smem, sdom->subDStream >>> 
	    ( nBlocks_j, nBlocks_k, i_inc, j_inc, k_inc, (kernel%(num_directions/FluxGroupDim))*FluxGroupDim, 
	      kernel/(num_directions/FluxGroupDim), num_groups, num_directions,
	      local_imax, local_jmax, local_kmax, d_dx, d_dy, d_dz, d_rhs, d_sigt, d_direction, d_psi, d_i_plane, d_j_plane, d_k_plane, kernel, nKernels);
  cudaCheckError();
	  }
	}
      }
      else {
	if ( k_inc == 0 ) {
	  if ( kernel%NumSMs == 0 ) {
  cudaCheckError();
	  sweep_over_hyperplane_ZGD_fluxRegisters <0,1,0>
	    <<< NumSMs, threadsPerBlock, required_smem, sdom->subDStream >>> 
	    ( nBlocks_j, nBlocks_k, i_inc, j_inc, k_inc, (kernel%(num_directions/FluxGroupDim))*FluxGroupDim, 
	      kernel/(num_directions/FluxGroupDim), num_groups, num_directions,
	      local_imax, local_jmax, local_kmax, d_dx, d_dy, d_dz, d_rhs, d_sigt, d_direction, d_psi, d_i_plane, d_j_plane, d_k_plane, kernel, nKernels);	
  cudaCheckError();
	  }
	}
	else {
	  if ( kernel%NumSMs == 0 ) {
  cudaCheckError();
	  sweep_over_hyperplane_ZGD_fluxRegisters <0,1,1>
	    <<< NumSMs, threadsPerBlock, required_smem, sdom->subDStream >>> 
	    ( nBlocks_j, nBlocks_k, i_inc, j_inc, k_inc, (kernel%(num_directions/FluxGroupDim))*FluxGroupDim, 
	      kernel/(num_directions/FluxGroupDim), num_groups, num_directions,
	      local_imax, local_jmax, local_kmax, d_dx, d_dy, d_dz, d_rhs, d_sigt, d_direction, d_psi, d_i_plane, d_j_plane, d_k_plane, kernel, nKernels);

	  cudaCheckError();

	  }
	}
      }
    }
    else {
      if ( j_inc == 0 ) {
	if ( k_inc == 0 ) {
	  if ( kernel%NumSMs == 0 ) {
  cudaCheckError();
	  sweep_over_hyperplane_ZGD_fluxRegisters <1,0,0>
	    <<< NumSMs, threadsPerBlock, required_smem, sdom->subDStream >>> 
	    ( nBlocks_j,
	      nBlocks_k,
	      i_inc,
	      j_inc,
	      k_inc,
	      (kernel%(num_directions/FluxGroupDim))*FluxGroupDim,
	      kernel/(num_directions/FluxGroupDim),
	      num_groups,
	      num_directions,
	      local_imax,
	      local_jmax,
	      local_kmax,
	      d_dx, 
	      d_dy, 
	      d_dz, 
	      d_rhs, 
	      d_sigt, 
	      d_direction,
	      d_psi, 
	      d_i_plane,                            
	      d_j_plane,                            
	      d_k_plane,
	      kernel,
	      nKernels
	      );
  cudaCheckError();
	  }
	}
	else {
	  if ( kernel%NumSMs == 0 ) {
  cudaCheckError();
	sweep_over_hyperplane_ZGD_fluxRegisters <1,0,1>
	  <<< NumSMs, threadsPerBlock, required_smem, sdom->subDStream >>> 
	  ( nBlocks_j, nBlocks_k, i_inc, j_inc, k_inc, (kernel%(num_directions/FluxGroupDim))*FluxGroupDim, 
	    kernel/(num_directions/FluxGroupDim), num_groups, num_directions,
	    local_imax, local_jmax, local_kmax, d_dx, d_dy, d_dz, d_rhs, d_sigt, d_direction, d_psi, d_i_plane, d_j_plane, d_k_plane, kernel, nKernels);	
  cudaCheckError();
	  }
	}
      }
      else {
	if ( k_inc == 0  ) {
	  if ( kernel%NumSMs == 0 ) {
  cudaCheckError();
	  sweep_over_hyperplane_ZGD_fluxRegisters <1,1,0>
	    <<< NumSMs, threadsPerBlock, required_smem, sdom->subDStream >>> 
	    ( nBlocks_j, nBlocks_k, i_inc, j_inc, k_inc, (kernel%(num_directions/FluxGroupDim))*FluxGroupDim, 
	      kernel/(num_directions/FluxGroupDim), num_groups, num_directions,
	      local_imax, local_jmax, local_kmax, d_dx, d_dy, d_dz, d_rhs, d_sigt, d_direction, d_psi, d_i_plane, d_j_plane, d_k_plane, kernel, nKernels);	
  cudaCheckError();
	  }
	}
	else {
	  if ( kernel%NumSMs == 0 ) {
  cudaCheckError();
	sweep_over_hyperplane_ZGD_fluxRegisters <1,1,1>
	  <<< NumSMs, threadsPerBlock, required_smem, sdom->subDStream >>> 
	  ( nBlocks_j, nBlocks_k, i_inc, j_inc, k_inc, (kernel%(num_directions/FluxGroupDim))*FluxGroupDim, 
	    kernel/(num_directions/FluxGroupDim), num_groups, num_directions,
	    local_imax, local_jmax, local_kmax, d_dx, d_dy, d_dz, d_rhs, d_sigt, d_direction, d_psi, d_i_plane, d_j_plane, d_k_plane, kernel, nKernels);	
  cudaCheckError();
	  }
	}
      }
    }
    
  }

  cudaCheckError();
  //cudaDeviceSynchronize();  
  omptime += omp_get_wtime();
  //  printf ("bandwidth achieved = %e GB/s, %d, %d, %d, %d, %d, %e \n",
  //	  16.0e-9/omptime*local_imax*local_jmax*local_kmax*num_directions*num_groups, local_imax, local_jmax, local_kmax, num_directions, num_groups, omptime);

  int ierr, my_id;

#ifdef OCTAVE_PROFILING
  ierr = MPI_Comm_rank (MPI_COMM_WORLD, &my_id);
  printf ("h = rectangle ('Position', [%f, %f, %f, %f]); set ( h, 'FaceColor', [%d, %d, %d]); \n", 
	  mpi_time_start, 1.0*my_id, MPI_Wtime()-mpi_time_start, 1.0, (i_inc==0) ? 0 : 1, (j_inc==0) ? 0 : 1, (k_inc==0) ? 0 : 1 );
#endif

  }

  /*
  printf ("\nKernel Time = %e s   ", time_s);
  printf ("Bandwidth achieved = %e GB/s \n", 
	  1.0e-9 * (                                                                               // data from GMEM
		    8.0 * local_imax * local_jmax * local_kmax * num_directions * num_groups       // rhs - read
		    + 8.0 * local_imax * local_jmax * local_kmax * num_directions * num_groups     // psi - write
		    + 8.0 * local_jmax*local_kmax * num_directions * num_groups * 2                // flux_i - r/w
		    + 8.0 * local_imax*local_kmax * num_directions * num_groups * 2 * local_jmax/8 // flux_j - r/w
		    + 8.0 * local_imax*local_jmax * num_directions * num_groups * 2 * local_kmax/8 // flux_k - r/w
		    )
	  / time_s);
  */
  
#ifdef CU_TIMING
  cudaEventRecord(start);
#endif

#ifdef CU_TIMING
  cudaEventRecord(stop);
  cudaDeviceSynchronize();
  cudaCheckError();
  cudaEventElapsedTime(&time_ms,start,stop);
  time_s=time_ms*.001;
  printf("ZGD: time to copy PSI D2H: %g [s]\n",time_s);
#endif

  return 0;

}


int cuda_sweep_ZGD( double *d_rhs, double *h_phi, double *h_psi, double *d_sigt, Directions *d_direction,
                    double *h_i_plane, double *h_j_plane, double *h_k_plane,
                    int *d_ii_jj_kk_z_idx, int *h_offset,  int *d_offset, double *d_dx, double *d_dy, double *d_dz,
                    int num_zones, int num_directions, int num_groups,
                    int local_imax, int local_jmax, int local_kmax, int Nslices){

  size_t N;
  size_t groups_dirs;
  double *d_phi;
  float time_ms, time_s;

  cudaCheckError();  

  N = num_zones * num_directions * num_groups;
  groups_dirs = num_directions * num_groups;


  cudaEvent_t start,stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);



#ifdef USE_PSI_HOST_MEM
  double *d_psi = h_psi;
#else
  double *d_psi;
  #ifdef CU_TIMING
  cudaEventRecord(start);
  #endif

  cudaMalloc((void **) &d_psi, N*sizeof(double));
  cudaMemcpy(d_psi, h_psi,   N*sizeof(double), cudaMemcpyHostToDevice);
  cudaCheckError();

  #ifdef CU_TIMING
  cudaEventRecord(stop);
  cudaDeviceSynchronize();
  cudaCheckError();
  cudaEventElapsedTime(&time_ms,start,stop);
  time_s=time_ms*.001;
  printf("ZGD: time to copy PSI H2D: %g [s]\n",time_s);
  #endif
#endif



  int i_plane_zones = local_jmax * local_kmax * groups_dirs;
  int j_plane_zones = local_imax * local_kmax * groups_dirs;
  int k_plane_zones = local_imax * local_jmax * groups_dirs;

#ifdef USE_IJK_PLANE_HOST_MEM

  double *d_i_plane = h_i_plane;
  double *d_j_plane = h_j_plane;
  double *d_k_plane = h_k_plane;

#else

  double *d_i_plane,  *d_j_plane, *d_k_plane;

  #ifdef CU_TIMING
  cudaEventRecord(start);
  #endif

  cudaMalloc((void **) &d_i_plane, i_plane_zones * sizeof(double));
  cudaMalloc((void **) &d_j_plane, j_plane_zones * sizeof(double));
  cudaMalloc((void **) &d_k_plane, k_plane_zones * sizeof(double));
  cudaMemcpy(d_i_plane, h_i_plane, i_plane_zones * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(d_j_plane, h_j_plane, j_plane_zones * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(d_k_plane, h_k_plane, k_plane_zones * sizeof(double), cudaMemcpyHostToDevice);
  cudaCheckError();

  #ifdef CU_TIMING
  cudaEventRecord(stop);
  cudaDeviceSynchronize();
  cudaCheckError();
  cudaEventElapsedTime(&time_ms,start,stop);
  time_s=time_ms*.001;
  printf("ZGD: time to copy ijk_PLANE H2D: %g [s]\n",time_s);
  #endif

#endif

  cudaCheckError();
 
//call cuda kernel to sweep over hyperplanes(slices)

  dim3 threadsPerBlock(32,8);

  for (int slice = 0; slice < Nslices; slice++){
    
     #ifdef CU_TIMING
     cudaEventRecord(start);                                             
     #endif

     dim3 numBlocks = h_offset[slice+1] - h_offset[slice];
     sweep_over_hyperplane_ZGD<<<numBlocks,threadsPerBlock>>>(slice,d_offset,d_ii_jj_kk_z_idx,num_directions,num_groups,local_imax,local_jmax,local_kmax,
                                                          d_dx, d_dy, d_dz, d_rhs, d_phi, d_psi, d_sigt, d_direction, d_i_plane, d_j_plane, d_k_plane);
     #ifdef CU_TIMING
     cudaEventRecord(stop);                                              
     cudaDeviceSynchronize();                                            
     cudaCheckError();                                                   
     float time_ms, time_s;                                              
     cudaEventElapsedTime(&time_ms,start,stop);                          
     time_s=time_ms*.001;                                                
     printf("ZGD: #blocks=%d, time=%g [s]\n",h_offset[slice+1] - h_offset[slice],time_s);
     #endif
     
  }

  #ifdef CU_TIMING
  cudaEventRecord(start);
  #endif

#ifndef USE_IJK_PLANE_HOST_MEM
    cudaMemcpy(h_i_plane, d_i_plane, i_plane_zones * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_j_plane, d_j_plane, j_plane_zones * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_k_plane, d_k_plane, k_plane_zones * sizeof(double), cudaMemcpyDeviceToHost);
#endif

  #ifdef CU_TIMING
  cudaEventRecord(stop);
  cudaDeviceSynchronize();
  cudaCheckError();
  cudaEventElapsedTime(&time_ms,start,stop);
  time_s=time_ms*.001;
  printf("ZGD: time to copy ijk_PLANE D2H: %g [s]\n",time_s);
  #endif


#ifndef USE_PSI_HOST_MEM

#ifdef CU_TIMING
  cudaEventRecord(start);
#endif

  cudaMemcpy(h_psi,     d_psi, N*sizeof(double),                   cudaMemcpyDeviceToHost);

#ifdef CU_TIMING
  cudaEventRecord(stop);
  cudaDeviceSynchronize();
  cudaCheckError();
  cudaEventElapsedTime(&time_ms,start,stop);
  time_s=time_ms*.001;
  printf("ZGD: time to copy PSI D2H: %g [s]\n",time_s);
#endif

  cudaFree(d_psi);

#endif


#ifndef USE_IJK_PLANE_HOST_MEM
  cudaFree(d_i_plane);
  cudaFree(d_j_plane);
  cudaFree(d_k_plane);
#endif

  cudaCheckError();

  return 0;
}


#define I_PLANE_INDEX(j, k) ((k)*(local_jmax) + (j))
#define J_PLANE_INDEX(i, k) ((k)*(local_imax) + (i))
#define K_PLANE_INDEX(i, j) ((j)*(local_imax) + (i))


__global__ void  sweep_over_hyperplane_ZGD(int sliceID, int * __restrict__ offset, int * __restrict__ ii_jj_kk_z_idx, int num_directions, int num_groups,
                    int local_imax, int local_jmax, int local_kmax,
                    double * __restrict__ dx, double * __restrict__ dy, double * __restrict__ dz,
                    double * __restrict__ rhs, double * __restrict__ phi, double * __restrict__ psi,
                    const double * __restrict__ sigt, const Directions * __restrict__ direction,
                    double *i_plane, double *j_plane, double *k_plane){

 
//each block will process one element 
      int element = offset[sliceID] + blockIdx.x;
      if (element > offset[sliceID+1]) return; 

      int i    = ii_jj_kk_z_idx[element*4];
      int j    = ii_jj_kk_z_idx[element*4+1];
      int k    = ii_jj_kk_z_idx[element*4+2];
      int z    = ii_jj_kk_z_idx[element*4+3];
      int dir_grp = num_directions*num_groups;
      int I_P_I = I_PLANE_INDEX(j, k);
      int J_P_I = J_PLANE_INDEX(i, k);
      int K_P_I = K_PLANE_INDEX(i, j);
      double two_inv_dxi = 2.0/dx[i + 1];
      double two_inv_dyj = 2.0/dy[j + 1];
      double two_inv_dzk = 2.0/dz[k + 1];
     

      double * KRESTRICT  block_rhs = &rhs[z*dir_grp];
//    double * KRESTRICT  block_phi = &phi[z*num_directions*num_groups];
      double * KRESTRICT  block_psi = &psi[z*dir_grp];
      const double * KRESTRICT  block_sigt = &sigt[z*num_groups];

      double * KRESTRICT psi_lf_z = &i_plane[I_P_I*dir_grp]; 
      double * KRESTRICT psi_fr_z = &j_plane[J_P_I*dir_grp]; 
      double * KRESTRICT psi_bo_z = &k_plane[K_P_I*dir_grp]; 

     

      for (int group = threadIdx.y; group < num_groups; group += blockDim.y){

          for (int  d = threadIdx.x; d < num_directions; d += blockDim.x){

            int gd = d + group*num_directions;

            double xcos_dxi =  direction[d].xcos * two_inv_dxi; 
            double ycos_dyj =  direction[d].ycos * two_inv_dyj;
            double zcos_dzk =  direction[d].zcos * two_inv_dzk;

            double psi_lf_z_g_d = psi_lf_z[gd];
            double psi_fr_z_g_d = psi_fr_z[gd];
            double psi_bo_z_g_d = psi_bo_z[gd];

            /* Calculate new zonal flux */
            double psi_z_g_d = (block_rhs[gd]
                + psi_lf_z_g_d * xcos_dxi
                + psi_fr_z_g_d * ycos_dyj
                + psi_bo_z_g_d * zcos_dzk)
                / (xcos_dxi + ycos_dyj + zcos_dzk + block_sigt[group]);

            block_psi[gd] = psi_z_g_d;


            /* Apply diamond-difference relationships */
            psi_lf_z[gd] = 2.0 * psi_z_g_d - psi_lf_z_g_d;
            psi_fr_z[gd] = 2.0 * psi_z_g_d - psi_fr_z_g_d;
            psi_bo_z[gd] = 2.0 * psi_z_g_d - psi_bo_z_g_d;
          }
        }
}

int cuda_LPlusTimes_sweep_ZGD( double *d_phi_out, double *d_ell_plus, double *h_psi, double *d_sigt, Directions *d_direction,
                    double *h_i_plane, double *h_j_plane, double *h_k_plane,
                    int *d_ii_jj_kk_z_idx, int *h_offset,  int *d_offset, double *d_dx, double *d_dy, double *d_dz,
                    int num_zones, int num_directions, int num_groups, int num_local_groups,
                    int nidx, int group0,
                    int local_imax, int local_jmax, int local_kmax, int Nslices){


  size_t N;
  size_t groups_dirs;
  float time_ms, time_s;
  
  N = num_zones * num_directions * num_local_groups;
  groups_dirs = num_directions * num_local_groups;

  cudaEvent_t start,stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);


#ifdef USE_PSI_HOST_MEM
  double *d_psi = h_psi;
#else
  double *d_psi;
  #ifdef CU_TIMING
  cudaEventRecord(start);
  #endif
  printf ("SHOULD NOT BE CALLING THIS \n");
  cudaMalloc((void **) &d_psi, N*sizeof(double));
  cudaMemcpy(d_psi, h_psi,   N*sizeof(double), cudaMemcpyHostToDevice);
  cudaCheckError();

  #ifdef CU_TIMING
  cudaEventRecord(stop);
  cudaDeviceSynchronize();
  cudaCheckError();
  cudaEventElapsedTime(&time_ms,start,stop);
  time_s=time_ms*.001;
  printf("ZGD: time to copy PSI H2D: %g [s]\n",time_s);
  #endif
#endif

  int i_plane_zones = local_jmax * local_kmax * groups_dirs;
  int j_plane_zones = local_imax * local_kmax * groups_dirs;
  int k_plane_zones = local_imax * local_jmax * groups_dirs;

#ifdef USE_IJK_PLANE_HOST_MEM

  double *d_i_plane = h_i_plane;
  double *d_j_plane = h_j_plane;
  double *d_k_plane = h_k_plane;

#else

  double *d_i_plane,  *d_j_plane, *d_k_plane;

  #ifdef CU_TIMING
  cudaEventRecord(start);
  #endif

  cudaMalloc((void **) &d_i_plane, i_plane_zones * sizeof(double));
  cudaMalloc((void **) &d_j_plane, j_plane_zones * sizeof(double));
  cudaMalloc((void **) &d_k_plane, k_plane_zones * sizeof(double));
  cudaMemcpy(d_i_plane, h_i_plane, i_plane_zones * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(d_j_plane, h_j_plane, j_plane_zones * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(d_k_plane, h_k_plane, k_plane_zones * sizeof(double), cudaMemcpyHostToDevice);
  cudaCheckError();

  #ifdef CU_TIMING
  cudaEventRecord(stop);
  cudaDeviceSynchronize();
  cudaCheckError();
  cudaEventElapsedTime(&time_ms,start,stop);
  time_s=time_ms*.001;
  printf("ZGD: time to copy ijk_PLANE H2D: %g [s]\n",time_s);
  #endif

#endif

  cudaCheckError();
 

  cudaFuncSetCacheConfig(sweep_over_hyperplane_ZGD, cudaFuncCachePreferL1);

//call cuda kernel to sweep over hyperplanes(slices)
  int dim_y = 8; 
  dim3 threadsPerBlock(32,dim_y);

  for (int slice = 0; slice < Nslices; slice++){
    
     #ifdef CU_TIMING
     cudaEventRecord(start);                                             
     #endif

     dim3 numBlocks = h_offset[slice+1] - h_offset[slice];
     LPlusTimes_sweep_over_hyperplane_ZGD<<<numBlocks,threadsPerBlock,num_directions*dim_y*sizeof(double)>>>(slice,d_offset,d_ii_jj_kk_z_idx,num_directions,num_groups,num_local_groups,nidx,group0, local_imax,local_jmax,local_kmax,
                                                          d_dx, d_dy, d_dz, d_phi_out, d_ell_plus, d_psi, d_sigt, d_direction, d_i_plane, d_j_plane, d_k_plane);
     #ifdef CU_TIMING
     cudaEventRecord(stop);                                              
     cudaDeviceSynchronize();                                            
     cudaCheckError();                                                   
     float time_ms, time_s;                                              
     cudaEventElapsedTime(&time_ms,start,stop);                          
     time_s=time_ms*.001;                                                
     printf("ZGD: #blocks=%d, time=%g [s]\n",h_offset[slice+1] - h_offset[slice],time_s);
     #endif
     
  }

  #ifdef CU_TIMING
  cudaEventRecord(start);
  #endif

#ifndef USE_IJK_PLANE_HOST_MEM
    cudaMemcpy(h_i_plane, d_i_plane, i_plane_zones * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_j_plane, d_j_plane, j_plane_zones * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_k_plane, d_k_plane, k_plane_zones * sizeof(double), cudaMemcpyDeviceToHost);
#endif

  #ifdef CU_TIMING
  cudaEventRecord(stop);
  cudaDeviceSynchronize();
  cudaCheckError();
  cudaEventElapsedTime(&time_ms,start,stop);
  time_s=time_ms*.001;
  printf("ZGD: time to copy ijk_PLANE D2H: %g [s]\n",time_s);
  #endif


#ifndef USE_PSI_HOST_MEM

#ifdef CU_TIMING
  cudaEventRecord(start);
#endif

  cudaMemcpy(h_psi,     d_psi, N*sizeof(double),                   cudaMemcpyDeviceToHost);

#ifdef CU_TIMING
  cudaEventRecord(stop);
  cudaDeviceSynchronize();
  cudaCheckError();
  cudaEventElapsedTime(&time_ms,start,stop);
  time_s=time_ms*.001;
  printf("ZGD: time to copy PSI D2H: %g [s]\n",time_s);
#endif

  cudaFree(d_psi);
#endif


//  cudaFree(d_phi);

#ifndef USE_IJK_PLANE_HOST_MEM
  cudaFree(d_i_plane);
  cudaFree(d_j_plane);
  cudaFree(d_k_plane);
#endif

  cudaCheckError();


  return 0;
}


__global__ void LPlusTimes_sweep_over_hyperplane_ZGD(int sliceID, int * __restrict__ offset, int * __restrict__ ii_jj_kk_z_idx, int num_directions, int num_groups, int num_local_groups, int nidx, int group0,
                    int local_imax, int local_jmax, int local_kmax, double * __restrict__ dx, double * __restrict__ dy, 
                    double * __restrict__ dz, double *__restrict__ phi_out, double * __restrict__ ell_plus, double * __restrict__ psi, 
                    const double * __restrict__ sigt, const Directions * __restrict__ direction,
                    double *i_plane, double *j_plane, double *k_plane){
 
     extern __shared__ double rhs_local[];
 
//each block will process one element 
      int element = offset[sliceID] + blockIdx.x;
      if (element > offset[sliceID+1]) return; 

      int i    = ii_jj_kk_z_idx[element*4];
      int j    = ii_jj_kk_z_idx[element*4+1];
      int k    = ii_jj_kk_z_idx[element*4+2];
      int z    = ii_jj_kk_z_idx[element*4+3];
      int dir_grp = num_directions*num_local_groups;
      int I_P_I = I_PLANE_INDEX(j, k);
      int J_P_I = J_PLANE_INDEX(i, k);
      int K_P_I = K_PLANE_INDEX(i, j);
      double two_inv_dxi = 2.0/dx[i + 1];
      double two_inv_dyj = 2.0/dy[j + 1];
      double two_inv_dzk = 2.0/dz[k + 1];
     
      double * KRESTRICT  block_psi = &psi[z*dir_grp];
      const double * KRESTRICT  block_sigt = &sigt[z*num_local_groups];
      double * KRESTRICT psi_lf_z = &i_plane[I_P_I*dir_grp]; 
      double * KRESTRICT psi_fr_z = &j_plane[J_P_I*dir_grp]; 
      double * KRESTRICT psi_bo_z = &k_plane[K_P_I*dir_grp];
 
      const double * KRESTRICT block_phi_out = &phi_out[z*num_groups*nidx + group0*nidx];

      for (int group = threadIdx.y; group < num_local_groups; group += blockDim.y){

	  double * __restrict__ rhs_local_group = &rhs_local[threadIdx.y * num_directions];
          for (int d = threadIdx.x; d < num_directions; d+=blockDim.x) {
            double sum  = 0.0;
            for(int nm_offset = 0;nm_offset < nidx;++nm_offset)
               sum += ell_plus[nm_offset + d*nidx] * block_phi_out[nm_offset + group*nidx];
            rhs_local_group[d] = sum;     
          }	  
	 
          for (int  d = threadIdx.x; d < num_directions; d += blockDim.x){

            int gd = d + group*num_directions;

            double xcos_dxi =  direction[d].xcos * two_inv_dxi; 
            double ycos_dyj =  direction[d].ycos * two_inv_dyj;
            double zcos_dzk =  direction[d].zcos * two_inv_dzk;

            double psi_lf_z_g_d = psi_lf_z[gd];
            double psi_fr_z_g_d = psi_fr_z[gd];
            double psi_bo_z_g_d = psi_bo_z[gd];

            /* Calculate new zonal flux */ 
            double psi_z_g_d = (rhs_local_group[d]
                + psi_lf_z_g_d * xcos_dxi
                + psi_fr_z_g_d * ycos_dyj 
                + psi_bo_z_g_d * zcos_dzk)
                / (xcos_dxi + ycos_dyj + zcos_dzk + block_sigt[group]);

            block_psi[gd] = psi_z_g_d;


            /* Apply diamond-difference relationships */
            psi_lf_z[gd] = 2.0 * psi_z_g_d - psi_lf_z_g_d;
            psi_fr_z[gd] = 2.0 * psi_z_g_d - psi_fr_z_g_d;
            psi_bo_z[gd] = 2.0 * psi_z_g_d - psi_bo_z_g_d;
          }
        }
}
