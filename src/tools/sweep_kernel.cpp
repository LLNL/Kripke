/*--------------------------------------------------------------------------
 * Main program for 3-D Sweep Kernel
 *--------------------------------------------------------------------------*/

#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#ifdef HAVE_PERFTOOLS_PROFILER
#include "Kripke/profiler.h"
#endif

#include <Kripke/transport_protos.h>
#include <Kripke/user_data.h>
#include <Kripke/comm.h>

int main(int argc, char *argv[])
{
  FILE             *in_file;
  FILE             *out_file;
  int myid;
  int ierr=0;
  int i;

  /*-----------------------------------------------------------------------
   * Initialize MPI
   *-----------------------------------------------------------------------*/
  MPI_Init(&argc, &argv);

  /*-----------------------------------------------------------------------
   * Open input and  output files
   *-----------------------------------------------------------------------*/
  MPI_Comm_rank(MPI_COMM_WORLD, &myid );

  if(myid == 0){
    /* Print out a banner message along with a version number. */
    printf("\n");
    printf("---------------------------------------------------------\n");
    printf("--------------- SWEEP KERNEL VERSION 1.0 ----------------\n");
    printf("---------------------------------------------------------\n");
  }

  if(argc < 2){
    printf("ERROR: Please specify an input file name.\n");
    error_exit(1);
  }
  std::string input_file_name = argv[1];


  /*-----------------------------------------------------------------------
   * Initialize user input and put relevant data in the user_data structure
   *-----------------------------------------------------------------------*/
  Input_Variables input_variables;
  input_variables.read(input_file_name);
  input_variables.print();
  User_Data *user_data = new User_Data(MPI_COMM_WORLD, &input_variables);

  /* Begin timing of solve */
  user_data->timing.start("Solve");

  /* Run the driver */
  for(i=0; i<user_data->ncalls; i++){
    SweepDriver(user_data);
  }

  /* End timing of solve */
  user_data->timing.stop("Solve");

  /*-----------------------------------------------------------------------
   * Print timing information
   *-----------------------------------------------------------------------*/
  user_data->timing.print();

  /*-----------------------------------------------------------------------
   * Cleanup and exist
   *-----------------------------------------------------------------------*/
  delete user_data;
  MPI_Finalize();

  return(0);
}
