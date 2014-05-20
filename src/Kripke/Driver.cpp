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


/*--------------------------------------------------------------------------
 * SweepDriver : Solves the sweeping problem defined in user_data
 *--------------------------------------------------------------------------*/

void Driver(User_Data *user_data)
{
  double sum=0.;
  double gsum=0.;
  double global_num_zones = user_data->global_num_zones;

  int *nzones = user_data->grid_data->nzones;
  int num_directions = user_data->directions.size();
  double total_groups = user_data->num_group_sets * user_data->num_groups_per_set;
  int num_zones = nzones[0]*nzones[1]*nzones[2];
  int d, zone;
  int myid;

  /* Begin timing of solve */
  user_data->timing.start("Solve");
  for(int iter = 0;iter < user_data->niter;++ iter){
    SweepSolver(user_data);
  }

  /* End timing of solve */
  user_data->timing.stop("Solve");
}
