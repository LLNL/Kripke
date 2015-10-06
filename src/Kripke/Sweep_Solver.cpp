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



//LG
#include<Kripke/SubTVec.h>

#ifdef KRIPKE_USE_CUDA
#include "Kripke/cu_utils.h"
#endif


#define KRIPKE_USE_CUDA_TIMING


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
      //MPI_Barrier(MPI_COMM_WORLD);
      BLOCK_TIMER(grid_data->timing, Source);
      kernel->source(grid_data);
      SYNC_GPU();
      //MPI_Barrier(MPI_COMM_WORLD);
    }

    // Moments to Discrete transformation
    {
      //MPI_Barrier(MPI_COMM_WORLD);
      BLOCK_TIMER(grid_data->timing, LPlusTimes);
      kernel->LPlusTimes(grid_data);
      SYNC_GPU();
      //MPI_Barrier(MPI_COMM_WORLD);
    }

    //grid_data->particleEdit();
    /*
     * Sweep each Group Set
     */
    { 
      //MPI_Barrier(MPI_COMM_WORLD);
      BLOCK_TIMER(grid_data->timing, Sweep);

#ifdef LPlusTimes_sweep_LTimes_combined
      #ifdef KRIPKE_USE_CUDA
      for(int ds = 0;ds < grid_data->num_zone_sets;++ ds)
        set_cudaMemZeroAsync( (void *) grid_data->d_phi[ds], (size_t)(grid_data->phi[ds]->elements) * sizeof(double));
      #endif
#endif
      if(true){
//      if(else){
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
      SYNC_GPU();
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
  // Create a new sweep communicator object
  SweepComm sweep_comm(grid_data);

  // Add all subdomains in our list
  for(int i = 0;i < subdomain_list.size();++ i){

    int sdom_id = subdomain_list[i];
    Subdomain &sdom = grid_data->subdomains[sdom_id];


    sweep_comm.addSubdomain(sdom_id, sdom);
    // Clear boundary conditions
    for(int dim = 0;dim < 3;++ dim){
      if(sdom.upwind[dim].subdomain_id == -1){
        sdom.plane_data[dim]->clear(0.0);
      }
    }
  }

  /* Loop until we have finished all of our work */
  while(sweep_comm.workRemaining()){

    // Get a list of subdomains that have met dependencies
    std::vector<int> sdom_ready = sweep_comm.readySubdomains();

    for(int idx = 0;idx < sdom_ready.size();++ idx){
      int sdom_id = sdom_ready[idx];

      /* Use standard Diamond-Difference sweep */
      {
        BLOCK_TIMER(grid_data->timing, Sweep_Kernel);

        Subdomain &sdom = grid_data->subdomains[sdom_id];
        grid_data->kernel->sweep(&sdom);
      }

      // Mark as complete (and do any communication)
      sweep_comm.markComplete(sdom_id);
    }
  }

  return(0);
}


