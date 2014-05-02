/*--------------------------------------------------------------------------
 * Driver routine for Sweep Kernel
 *--------------------------------------------------------------------------*/

#include "transport_headers.h"

/*--------------------------------------------------------------------------
 * SweepDriver : Solves the sweeping problem defined in user_data
 *--------------------------------------------------------------------------*/

void SweepDriver(User_Data *user_data)
{
  double **rhs, **psi, **tempv;
  double sum=0.;
  double gsum=0.;
  double global_num_zones = user_data->grid_data->global_num_zones;

  int *nzones = user_data->grid_data->nzones;
  int num_directions = user_data->grid_data->directions.size();
  int num_zones = nzones[0]*nzones[1]*nzones[2];
  int d, zone;
  int myid;

  /* Allocate space for rhs, psi, tempv */
  CNEW(rhs, num_directions, double **);
  CNEW(psi, num_directions, double **);
  CNEW(tempv, num_directions, double **);
  for(d=0; d<num_directions; d++){
    NEW(rhs[d], num_zones, double *);
    for(zone=0; zone<num_zones; zone++){
      rhs[d][zone] = user_data->source_value;
    }
    CNEW(psi[d], num_zones, double *);
    CNEW(tempv[d], num_zones, double *);
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
    FREE(rhs[d]);
    FREE(psi[d]);
    FREE(tempv[d]);
  }
  FREE(rhs);
  FREE(psi);
  FREE(tempv);

}
