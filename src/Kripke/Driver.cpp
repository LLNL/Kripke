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

  myid = GetRrank();
/*  if(myid == 0){
    // Print out a banner message along with a version number.
    printf("\n");
    printf("---------------------------------------------------------\n");
    printf("Performing Sweeps\n");
    printf("---------------------------------------------------------\n");
  }
*/
  for(int iter = 0;iter < user_data->niter;++ iter){
    if(myid == 0){
      printf("  iter %3d\n", iter);
    }
    SweepSolver(user_data);
  }

  /* Sum all entries in psi and output average */
  /*
  sum = 0.0;
  Grid_Data *grid_data = user_data->grid_data;
  for(int gs = 0;gs < grid_data->gd_sets.size();++ gs){
    for(int ds = 0;ds < grid_data->gd_sets[gs].size();++ ds){
      sum += grid_data->gd_sets[gs][ds].psi->sum();
    }
  }
  gsum = sum;
  MPI_Allreduce( &sum, &gsum, 1, MPI_DOUBLE, MPI_SUM, GetRGroup());
  sum = gsum/(global_num_zones*((double)num_directions)*total_groups);
  if(myid == 0){
    printf("\n");
    printf("Global number of zones = %22.16e\n",global_num_zones);
    printf("Global number of directions = %d\n", num_directions);
    printf("Global number of groups = %d\n", (int)total_groups);
    printf("Global sum = %22.16e\n", gsum);
    printf("Global sum ratio (should equal 1) = %22.16e\n", sum);
  }*/

  /* End timing of solve */
  user_data->timing.stop("Solve");
}
