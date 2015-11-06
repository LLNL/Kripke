/*--------------------------------------------------------------------------
 * Sweep-based solver routine.
 *--------------------------------------------------------------------------*/

#include <Kripke.h>
#include <Kripke/Subdomain.h>
#include <Kripke/SubTVec.h>
#include <Kripke/SweepComm.h>
#include <Kripke/Grid.h>
#include <vector>
#include <stdio.h>

#include <cuda.h>
#include <cuda_runtime.h>
#include "nvToolsExt.h"


//LG
#include<Kripke/SubTVec.h>

#ifdef KRIPKE_USE_CUDA
#include "Kripke/cu_utils.h"
#endif

cudaStream_t GlobalStreams[34];

//#define KRIPKE_USE_CUDA_TIMING

//Macro for checking cuda errors following a cuda launch or api call
#define cudaCheckError() {                                              \
    cudaError_t e=cudaGetLastError();					\
    if(e!=cudaSuccess) {						\
      printf("Cuda failure %s:%d: '%s'\n",__FILE__,__LINE__,cudaGetErrorString(e)); \
      exit(EXIT_FAILURE);						\
    }									\
  }

void SYNC_GPU(){
#ifdef KRIPKE_USE_CUDA

#ifdef KRIPKE_USE_CUDA_TIMING
  do_cudaDeviceSynchronize();
#endif
  ; 
#else
  ;
#endif
}


#include <stddef.h>
#include <sys/time.h>
const uint32_t colors[] = { 0x0000ff00, 0x000000ff, 0x00ffff00, 0x00ff00ff, 0x0000ffff, 0x00ff0000, 0x00ffffff };
const int num_colors = sizeof(colors)/sizeof(uint32_t);

#define PUSH_RANGE(name,cid) {				\
    int color_id = cid;					\
    color_id = color_id%num_colors;			\
    nvtxEventAttributes_t eventAttrib = {0};		\
    eventAttrib.version = NVTX_VERSION;			\
    eventAttrib.size = NVTX_EVENT_ATTRIB_STRUCT_SIZE;	\
    eventAttrib.colorType = NVTX_COLOR_ARGB;		\
    eventAttrib.color = colors[color_id];		\
    eventAttrib.messageType = NVTX_MESSAGE_TYPE_ASCII;	\
    eventAttrib.message.ascii = name;			\
    nvtxRangePushEx(&eventAttrib);			\
  }


template <int iOct, int jOct, int kOct>
__global__ void sweep_over_hyperplane_ZGD_fluxRegisters ( const int nBlocks_j,
							  const int nBlocks_k,
							  const int i_inc,
							  const int j_inc,
							  const int k_inc,
							  const int direction_offset,
							  const int group, 
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
							  double * __restrict__ flux_boundary_k
							  );

double mysecond (void)
{
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return (double)tv.tv_sec + (double)tv.tv_usec * 1.0e-6;
}


/**
   Run solver iterations.
*/
int SweepSolver (Grid_Data *grid_data)
{
  Kernel *kernel = grid_data->kernel;

  int mpi_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

  BLOCK_TIMER(grid_data->timing, Solve);

  for ( int i=0; i<34; i++ ) {
    cudaStreamCreate ( &(GlobalStreams[i]) );
    cudaCheckError();
  }
  // Loop over iterations
  double part_last = 0.0;
  for(int iter = 0;iter < grid_data->niter;++ iter){

    /*
     * Compute the RHS
     */

    // Discrete to Moments transformation
   
    {
      BLOCK_TIMER(grid_data->timing, LTimes);
#ifdef LPlusTimes_sweep_LTimes_combined
      if (iter == 0)
#endif
        kernel->LTimes(grid_data);
      SYNC_GPU();
    }

    // Compute Scattering Source Term
    {
      BLOCK_TIMER(grid_data->timing, Scattering);
      kernel->scattering(grid_data);
      SYNC_GPU();
    }

    // Compute External Source Term
    {
      BLOCK_TIMER(grid_data->timing, Source);
      kernel->source(grid_data);
      SYNC_GPU();
    }

    // Moments to Discrete transformation
    {
      BLOCK_TIMER(grid_data->timing, LPlusTimes);
      kernel->LPlusTimes(grid_data);
      SYNC_GPU();
    }

    //grid_data->particleEdit();
    /*
     * Sweep each Group Set
     */
    { 
      BLOCK_TIMER(grid_data->timing, Sweep);

#ifdef LPlusTimes_sweep_LTimes_combined
#ifdef KRIPKE_USE_CUDA
      for(int ds = 0;ds < grid_data->num_zone_sets;++ ds)
        set_cudaMemZeroAsync( (void *) grid_data->d_phi[ds], (size_t)(grid_data->phi[ds]->elements) * sizeof(double));
#endif
#endif
      if(true){
        // Create a list of all groups
        std::vector<int> sdom_list(grid_data->subdomains.size());
        for(int i = 0;i < grid_data->subdomains.size();++ i){
          sdom_list[i] = i;
        }

        // Sweep everything
        SweepSubdomains(sdom_list, grid_data);
      }
      // This is the ARDRA version, doing each groupset sweep independently
      else{
        for(int group_set = 0;group_set < grid_data->num_group_sets;++ group_set){
          std::vector<int> sdom_list;
          // Add all subdomains for this groupset
          for(int s = 0;s < grid_data->subdomains.size();++ s){
            if(grid_data->subdomains[s].idx_group_set == group_set){
              sdom_list.push_back(s);
            }
          }

          // Sweep the groupset
          SweepSubdomains(sdom_list, grid_data);
        }
      }
#ifdef KRIPKE_USE_CUDA_TIMING
      MPI_Barrier(MPI_COMM_WORLD);
      //SYNC_GPU();
#endif
    }

    {
#ifdef KRIPKE_USE_CUDA_TIMING
      MPI_Barrier(MPI_COMM_WORLD);
#endif
      BLOCK_TIMER(grid_data->timing, particleEdit);

      double part = grid_data->particleEdit();
      if(mpi_rank==0){
        printf("iter %d: particle count=%e, change=%e\n", iter, part, (part-part_last)/part);
      }

      part_last = part;
#ifdef KRIPKE_USE_CUDA_TIMING
      MPI_Barrier(MPI_COMM_WORLD);
#endif

    }

  }
  return(0);
}



/**
   Perform full parallel sweep algorithm on subset of subdomains.
*/
int SweepSubdomains (std::vector<int> subdomain_list, Grid_Data *grid_data)
{
  //printf ("Entering Sweep Subdomains \n");
  // Create a new sweep communicator object
  SweepComm sweep_comm(grid_data);

  // Add all subdomains in our list
  for(int i = 0;i < subdomain_list.size();++ i){

    int sdom_id = subdomain_list[i];
    Subdomain &sdomA = grid_data->subdomains[sdom_id];

#ifdef KRIPKE_USE_CUDA

    if(grid_data->kernel->sweep_mode == SWEEP_GPU){
      //LG  copy RHS to device
      double *dptr_h_rhs = sdomA.rhs->ptr();
      if ( sdomA.d_rhs == NULL){ // allocate
	sdomA.d_rhs = (double *) get_cudaMalloc((size_t) ( sdomA.num_zones * sdomA.num_groups * sdomA.num_directions) * sizeof(double));
      }
    }

#ifdef KRIPKE_ZGD_FLUX_REGISTERS

    // SR - copy up all required data before entering sweep loop
    if(grid_data->kernel->sweep_mode == SWEEP_GPU){

      cudaCheckError();

      Subdomain * sdom = &sdomA;

      int num_directions = sdom->num_directions;
      int num_groups = sdom->num_groups;
      int num_zones = sdom->num_zones;
      
      Directions *direction = sdom->directions;

      int local_imax = sdom->nzones[0];
      int local_jmax = sdom->nzones[1];
      int local_kmax = sdom->nzones[2];
      
      if ( sdom->d_psi == NULL ) {
	sdom->d_psi = (double *) get_cudaMalloc ( (size_t) (num_zones * num_directions * num_groups )*sizeof(double) );
	cudaCheckError();
      }
      
      if ( sdom->d_i_plane == NULL ) {

	cudaCheckError();

	size_t groups_dirs = num_directions * num_groups;
	sdom->d_i_plane = (double *) get_cudaMalloc ( (size_t) (local_jmax * local_kmax * groups_dirs )*sizeof(double) );
	sdom->d_j_plane = (double *) get_cudaMalloc ( (size_t) (local_imax * local_kmax * groups_dirs )*sizeof(double) );
	sdom->d_k_plane = (double *) get_cudaMalloc ( (size_t) (local_imax * local_jmax * groups_dirs )*sizeof(double) );
	cudaMemset (sdom->d_i_plane, 0, local_jmax * local_kmax * groups_dirs * sizeof(double) );
	cudaCheckError();
	cudaMemset (sdom->d_j_plane, 0, local_imax * local_kmax * groups_dirs * sizeof(double) );
	cudaCheckError();
	cudaMemset (sdom->d_k_plane, 0, local_imax * local_jmax * groups_dirs * sizeof(double) );

	cudaCheckError();

      }

      if ( sdom->d_streams == NULL ) {
	for ( int i=0; i<34; i++ ) {
	  sdom->sweepStreams[i] = GlobalStreams[i]; 
	}
	sdom->d_streams = new (double);
      }

      if ( sdom->d_events == NULL ) {
	for ( int i=0; i<34; i++ ) {
	  cudaEventCreate(&(sdom->sweepEvents[i]));
	  cudaCheckError();
	}
	sdom->d_events = new (double);
      }

      if ( sdom->tempi == NULL ) {
	size_t groups_dirs = num_directions * num_groups;
	int i_plane_zones = local_jmax * local_kmax * groups_dirs;
	int j_plane_zones = local_imax * local_kmax * groups_dirs;
	int k_plane_zones = local_imax * local_jmax * groups_dirs; 
	cudaMallocHost ( &(sdom->tempi), i_plane_zones * sizeof(double));
	cudaMallocHost ( &(sdom->tempj), j_plane_zones * sizeof(double));
	cudaMallocHost ( &(sdom->tempk), k_plane_zones * sizeof(double));
      }

      double *d_psi = sdom->d_psi;
      size_t N = num_zones * num_directions * num_groups;

      size_t groups_dirs = num_directions * num_groups;
      int i_plane_zones = local_jmax * local_kmax * groups_dirs;
      int j_plane_zones = local_imax * local_kmax * groups_dirs;
      int k_plane_zones = local_imax * local_jmax * groups_dirs; 
      
      double *d_i_plane,  *d_j_plane, *d_k_plane;

      d_i_plane = sdom->d_i_plane;
      d_j_plane = sdom->d_j_plane;
      d_k_plane = sdom->d_k_plane;

      double * tempi = sdom->tempi;
      double * tempj = sdom->tempj; 
      double * tempk = sdom->tempk;

      double *h_i_plane = sdom->plane_data[0]->ptr();
      double *h_j_plane = sdom->plane_data[1]->ptr();
      double *h_k_plane = sdom->plane_data[2]->ptr();

      cudaCheckError();

      cudaFuncSetCacheConfig ( sweep_over_hyperplane_ZGD_fluxRegisters<0,0,0>, cudaFuncCachePreferShared );
      cudaFuncSetCacheConfig ( sweep_over_hyperplane_ZGD_fluxRegisters<0,0,1>, cudaFuncCachePreferShared );
      cudaFuncSetCacheConfig ( sweep_over_hyperplane_ZGD_fluxRegisters<1,0,0>, cudaFuncCachePreferShared );
      cudaFuncSetCacheConfig ( sweep_over_hyperplane_ZGD_fluxRegisters<1,0,1>, cudaFuncCachePreferShared );
      cudaFuncSetCacheConfig ( sweep_over_hyperplane_ZGD_fluxRegisters<0,1,0>, cudaFuncCachePreferShared );
      cudaFuncSetCacheConfig ( sweep_over_hyperplane_ZGD_fluxRegisters<0,1,1>, cudaFuncCachePreferShared );
      cudaFuncSetCacheConfig ( sweep_over_hyperplane_ZGD_fluxRegisters<1,1,0>, cudaFuncCachePreferShared );
      cudaFuncSetCacheConfig ( sweep_over_hyperplane_ZGD_fluxRegisters<1,1,1>, cudaFuncCachePreferShared );

    }
#endif

#else

#endif
    
    sweep_comm.addSubdomain(sdom_id, sdomA);
    // Clear boundary conditions
    for(int dim = 0;dim < 3;++ dim){
      if(sdomA.upwind[dim].subdomain_id == -1){
        sdomA.plane_data[dim]->clear(0.0);
      }
    }
  }

  /* Loop until we have finished all of our work */
  while(sweep_comm.workRemaining()){

#ifdef KRIPKE_ZGD_FLUX_REGISTERS
   
    int  ierr, my_id;

    if(grid_data->kernel->sweep_mode == SWEEP_GPU){
      cudaCheckError();
      nvtxRangePushA ("readySubdomains");
    }

#endif

    std::vector<int> sdom_ready = sweep_comm.readySubdomains();

#ifdef KRIPKE_ZGD_FLUX_REGISTERS

    if(grid_data->kernel->sweep_mode == SWEEP_GPU){

      nvtxRangePop();

      cudaCheckError();

      // copy up first boundary planes
      if ( sdom_ready.size() > 0 ) {

	// ADD boundary plane data here
	int sdom_id = sdom_ready[0];
	Subdomain &sdomA = grid_data->subdomains[sdom_id];
	Subdomain * sdom = &sdomA;
	int num_directions = sdom->num_directions;
	int num_groups = sdom->num_groups;
	int local_imax = sdom->nzones[0];
	int local_jmax = sdom->nzones[1];
	int local_kmax = sdom->nzones[2];
	double *h_i_plane = sdom->plane_data[0]->ptr();
	double *h_j_plane = sdom->plane_data[1]->ptr();
	double *h_k_plane = sdom->plane_data[2]->ptr();
	double *d_i_plane,  *d_j_plane, *d_k_plane;
	d_i_plane = sdom->d_i_plane;
	d_j_plane = sdom->d_j_plane;
	d_k_plane = sdom->d_k_plane;
	size_t groups_dirs = num_directions * num_groups;
	int i_plane_zones = local_jmax * local_kmax * groups_dirs;
	int j_plane_zones = local_imax * local_kmax * groups_dirs;
	int k_plane_zones = local_imax * local_jmax * groups_dirs; 
	cudaCheckError();
	cudaMemcpyAsync(d_i_plane, h_i_plane, i_plane_zones * sizeof(double), cudaMemcpyHostToDevice, sdom->sweepStreams[32]);
	cudaCheckError();
	cudaMemcpyAsync(d_j_plane, h_j_plane, j_plane_zones * sizeof(double), cudaMemcpyHostToDevice, sdom->sweepStreams[32]);
	cudaCheckError();
	cudaMemcpyAsync(d_k_plane, h_k_plane, k_plane_zones * sizeof(double), cudaMemcpyHostToDevice, sdom->sweepStreams[32]);
	cudaCheckError();

      }

      cudaDeviceSynchronize();

      nvtxRangePushA("sweepLoop");
      cudaCheckError();

    }
#endif

    for(int idx = 0;idx < sdom_ready.size();++ idx){
      double mpi_time_start = MPI_Wtime();
      int sdom_id = sdom_ready[idx];

#ifdef KRIPKE_ZGD_FLUX_REGISTERS

      Subdomain * sdom;
      if(grid_data->kernel->sweep_mode == SWEEP_GPU){
	cudaCheckError();

	Subdomain &sdomA = grid_data->subdomains[sdom_id];
	sdom = &sdomA;
	for ( int i=0; i<34; i++ ) {
	  cudaStreamWaitEvent( sdom->sweepStreams[i], sdom->sweepEvents[0], 0 );
	}

	// copy up boundary data for NEXT subdomain - overlaps with compute
	if ( idx+1 < sdom_ready.size() ) {
	  int sdom_id = sdom_ready[idx+1];
	  Subdomain &sdomA = grid_data->subdomains[sdom_id];
	  Subdomain * sdom = &sdomA;
	  int num_directions = sdom->num_directions;
	  int num_groups = sdom->num_groups;
	  int local_imax = sdom->nzones[0];
	  int local_jmax = sdom->nzones[1];
	  int local_kmax = sdom->nzones[2];
	  double *h_i_plane = sdom->plane_data[0]->ptr();
	  double *h_j_plane = sdom->plane_data[1]->ptr();
	  double *h_k_plane = sdom->plane_data[2]->ptr();
	  double *d_i_plane,  *d_j_plane, *d_k_plane;
	  d_i_plane = sdom->d_i_plane;
	  d_j_plane = sdom->d_j_plane;
	  d_k_plane = sdom->d_k_plane;
	  size_t groups_dirs = num_directions * num_groups;
	  int i_plane_zones = local_jmax * local_kmax * groups_dirs;
	  int j_plane_zones = local_imax * local_kmax * groups_dirs;
	  int k_plane_zones = local_imax * local_jmax * groups_dirs; 
	  cudaMemcpyAsync(d_i_plane, h_i_plane, i_plane_zones * sizeof(double), cudaMemcpyHostToDevice, sdom->sweepStreams[32]);
	  cudaMemcpyAsync(d_j_plane, h_j_plane, j_plane_zones * sizeof(double), cudaMemcpyHostToDevice, sdom->sweepStreams[32]);
	  cudaMemcpyAsync(d_k_plane, h_k_plane, k_plane_zones * sizeof(double), cudaMemcpyHostToDevice, sdom->sweepStreams[32]);
	  cudaEventRecord( sdom->sweepEvents[0], sdom->sweepStreams[32] );
	}

#ifdef OCTAVE_PROFILING
	ierr = MPI_Comm_rank (MPI_COMM_WORLD, &my_id);
	printf ("h = rectangle ('Position', [%f, %f, %f, %f]); set ( h, 'FaceColor', [%f, %f, %f]); set ( h, 'EdgeColor', [%f, %f, %f]);  \n", 
		mpi_time_start, 1.0*my_id, MPI_Wtime()-mpi_time_start, 1.0, .9,.9,.9, 1.0, 0.0, 0.0 );
#endif

	cudaCheckError();
      }
#endif

      /* Use standard Diamond-Difference sweep */
      {
        BLOCK_TIMER(grid_data->timing, Sweep_Kernel);

        Subdomain &sdom = grid_data->subdomains[sdom_id];
        grid_data->kernel->sweep(&sdom);                  // SWEEP 
      }

#ifdef KRIPKE_ZGD_FLUX_REGISTERS

      if(grid_data->kernel->sweep_mode == SWEEP_GPU){

	if ( idx > 0 ) {
	  // So long as this is not the first domain, then a previous domain is nearing completion
	  // So wait for the D2H stream to clear and then mark the previous domain as complete
	  int sdom_id = sdom_ready[idx-1];
	  Subdomain &sdomA = grid_data->subdomains[sdom_id];
	  Subdomain * sdom = &sdomA;
	  cudaStreamSynchronize(sdom->sweepStreams[33]);
	  sweep_comm.markComplete(sdom_id);
	}

	cudaDeviceSynchronize();
      
	// copy down the boundary data.  no need to unswizzle.  asynchronously in dedicated copy down stream
	mpi_time_start = MPI_Wtime();
      
	int num_directions = sdom->num_directions;
	int num_groups = sdom->num_groups;
	int local_imax = sdom->nzones[0];
	int local_jmax = sdom->nzones[1];
	int local_kmax = sdom->nzones[2];
	double *h_i_plane = sdom->plane_data[0]->ptr();
	double *h_j_plane = sdom->plane_data[1]->ptr();
	double *h_k_plane = sdom->plane_data[2]->ptr();
	double *d_i_plane,  *d_j_plane, *d_k_plane;
	d_i_plane = sdom->d_i_plane;
	d_j_plane = sdom->d_j_plane;
	d_k_plane = sdom->d_k_plane;
	size_t groups_dirs = num_directions * num_groups;
	int i_plane_zones = local_jmax * local_kmax * groups_dirs;
	int j_plane_zones = local_imax * local_kmax * groups_dirs;
	int k_plane_zones = local_imax * local_jmax * groups_dirs; 
	cudaMemcpyAsync(h_i_plane, d_i_plane, i_plane_zones * sizeof(double), cudaMemcpyDeviceToHost, sdom->sweepStreams[33]);
	cudaMemcpyAsync(h_j_plane, d_j_plane, j_plane_zones * sizeof(double), cudaMemcpyDeviceToHost, sdom->sweepStreams[33]);
	cudaMemcpyAsync(h_k_plane, d_k_plane, k_plane_zones * sizeof(double), cudaMemcpyDeviceToHost, sdom->sweepStreams[33]);

#ifdef OCTAVE_PROFILING
	ierr = MPI_Comm_rank (MPI_COMM_WORLD, &my_id);
	printf ("h = rectangle ('Position', [%f, %f, %f, %f]); set ( h, 'FaceColor', [%f, %f, %f]); set ( h, 'EdgeColor', [%f, %f, %f]);  \n", 
		mpi_time_start, 1.0*my_id, MPI_Wtime()-mpi_time_start, 1.0, .8,.8,.8, 0.0, 1.0, 0.0 );
#endif

	// Mark as complete (and do any communication)
	if ( idx+1 ==  sdom_ready.size() ){
	  cudaStreamSynchronize(sdom->sweepStreams[33]);
	  sweep_comm.markComplete(sdom_id);
	}
      }
      else {

	// Mark as complete (and do any communication)
	sweep_comm.markComplete(sdom_id);

      }

#else

      // Mark as complete (and do any communication)
      sweep_comm.markComplete(sdom_id);

#endif 

    }

  }

#ifdef KRIPKE_ZGD_FLUX_REGISTERS


  if(grid_data->kernel->sweep_mode == SWEEP_GPU){
    nvtxRangePop();
      
    // Add all subdomains in our list
    for(int i = 0;i < subdomain_list.size();++ i){
      
      int sdom_id = subdomain_list[i];
      Subdomain &sdomA = grid_data->subdomains[sdom_id];
      
      // SR - copy back Psi data      
      Subdomain * sdom = &sdomA;
      cudaMemcpyAsync(sdom->psi->ptr(), sdom->d_psi, sdom->num_groups * sdom->num_zones * sdom->num_directions * sizeof(double), cudaMemcpyDeviceToHost, sdom->sweepStreams[0] );
      
    }
    
  }

#endif

  return(0);
}


