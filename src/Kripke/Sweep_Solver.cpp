/*--------------------------------------------------------------------------
 * Sweep-based solver routine.
 *--------------------------------------------------------------------------*/

#include <Kripke.h>
#include <Kripke/Comm.h>
#include <Kripke/Grid.h>
#include <vector>
#include <stdio.h>



/*----------------------------------------------------------------------
 * SweepSolverSolve
 *----------------------------------------------------------------------*/

int SweepSolver (Grid_Data *grid_data)
{
  Kernel *kernel = grid_data->kernel;

  BLOCK_TIMER(grid_data->timing, Solve);

  // Loop over iterations
  for(int iter = 0;iter < grid_data->niter;++ iter){

    /*
     * Compute the RHS
     */

    // Discrete to Moments transformation
    {
      BLOCK_TIMER(grid_data->timing, LTimes);
      kernel->LTimes(grid_data);
    }


    // This is where the Scattering kernel would go!



    // Moments to Discrete transformation
    {
      BLOCK_TIMER(grid_data->timing, LPlusTimes);
      kernel->LPlusTimes(grid_data);
    }

    /*
     * Sweep each Group Set
     */
    {
      BLOCK_TIMER(grid_data->timing, Sweep);
      for(int group_set = 0;group_set < grid_data->num_group_sets;++ group_set){
        SweepSolver_GroupSet(group_set, grid_data);
      }
    }
  }
  return(0);
}


/*----------------------------------------------------------------------
 * SweepSolverSolveDD
 *----------------------------------------------------------------------*/

int SweepSolver_GroupSet (int group_set, Grid_Data *grid_data)
{
  SweepComm sweep_comm;

  // Add all subdomains for this groupset
  for(int s = 0;s < grid_data->subdomains.size();++ s){
    if(grid_data->subdomains[s].idx_group_set == group_set){
      sweep_comm.addSubdomain(s, grid_data->subdomains[s]);
    }
  }

  /* Loop until we have finished all of our work */
  while(sweep_comm.workRemaining()){

    std::vector<int> sdom_ready = sweep_comm.readySubdomains();

    for(int idx = 0;idx < sdom_ready.size();++ idx){
      int sdom_id = sdom_ready[idx];

      /* Use standard Diamond-Difference sweep */
      {
        BLOCK_TIMER(grid_data->timing, Sweep_Kernel);

        Subdomain &sdom = grid_data->subdomains[sdom_id];
        grid_data->kernel->sweep(grid_data, &sdom, &sdom.plane_data[0][0], &sdom.plane_data[1][0], &sdom.plane_data[2][0]);
      }

      // Mark as complete (and do any communication)
      sweep_comm.markComplete(sdom_id);
    }
  }

  return(0);
}

