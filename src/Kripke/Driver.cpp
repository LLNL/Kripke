/*--------------------------------------------------------------------------
 * Driver routine for Sweep Kernel
 *--------------------------------------------------------------------------*/

#include <Kripke.h>
#include <Kripke/Comm.h>
#include <Kripke/Directions.h>
#include <Kripke/Grid.h>
#include <Kripke/User_Data.h>
#include <Kripke/SubTVec.h>
#include <stdio.h>


/**
 * Solves the sweeping problem defined in user_data.
 * Calls the SweepSolver niter times.
 */
void Driver(User_Data *user_data)
{
  BLOCK_TIMER(user_data->timing, Solve);

  // Loop over iterations
  for(int iter = 0;iter < user_data->niter;++ iter){
    SweepSolver(user_data);
  }
}
