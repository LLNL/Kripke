/*--------------------------------------------------------------------------
 * Driver routine for Sweep Kernel
 *--------------------------------------------------------------------------*/

#include <Kripke/transport_protos.h>
#include <Kripke/comm.h>
#include <Kripke/directions.h>
#include <Kripke/grid.h>
#include <Kripke/user_data.h>

#include<stdio.h>


/*--------------------------------------------------------------------------
 * SweepDriver : Solves the sweeping problem defined in user_data
 *--------------------------------------------------------------------------*/

void SweepDriver(User_Data *user_data)
{
  double sum=0.;
  double gsum=0.;
  double global_num_zones = user_data->global_num_zones;

  int *nzones = user_data->grid_data->nzones;
  int num_directions = user_data->directions.size();
  int num_zones = nzones[0]*nzones[1]*nzones[2];
  int d, zone;
  int myid;

  /* Allocate space for rhs, psi, tempv */
  double **rhs = new double*[num_directions];
  double **psi = new double*[num_directions];
  double **tempv = new double*[num_directions];
  for(d=0; d<num_directions; d++){
    rhs[d] = new double[num_zones];
    psi[d] = new double[num_zones];
    tempv[d] = new double[num_zones];

    for(zone=0; zone<num_zones; zone++){
      rhs[d][zone] = user_data->source_value;
      psi[d][zone] = 0.0;
      tempv[d][zone] = 0.0;
    }
  }

  myid = GetRrank();
  if(myid == 0){
    /* Print out a banner message along with a version number. */
    printf("\n");
    printf("---------------------------------------------------------\n");
    printf("Performing Sweeps\n");
    printf("---------------------------------------------------------\n");
  }

  /*---------------------------------------------------------------------
   * Call BoltzmannSolverSolve to solve the H Psi = R linear system
   * for all directions.
   *---------------------------------------------------------------------*/
  SweepSolverSolve(user_data, rhs, psi, tempv);

  /* Sum all entries in psi and output average */
  for(d=0; d<num_directions; d++){
    for(zone=0; zone<num_zones; zone++){
      sum += psi[d][zone];
    }
  }
  gsum = sum;
  MPI_Allreduce( &sum, &gsum, 1, MPI_DOUBLE, MPI_SUM, GetRGroup());
  sum = gsum/(global_num_zones*((double)num_directions));
  if(myid == 0){
    printf("\n");
    printf("Global number of zones = %22.16e\n",global_num_zones);
    printf("Global number of directions = %d\n", num_directions);
    printf("Global sum = %22.16e\n", gsum);
    printf("Global sum ratio (should equal 1) = %22.16e\n", sum);
  }

  /* Free rhs, psi, tempv */
  for(d=0; d<num_directions; d++){
    delete[] rhs[d];
    delete[] psi[d];
    delete[] tempv[d];
  }
  delete[] rhs;
  delete[] psi;
  delete[] tempv;

}
