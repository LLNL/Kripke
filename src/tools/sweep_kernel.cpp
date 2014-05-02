/*--------------------------------------------------------------------------
 * Main program for 3-D Sweep Kernel
 *--------------------------------------------------------------------------*/

#include<stdio.h>
#include<stdlib.h>

#ifdef HAVE_PERFTOOLS_PROFILER
#include "Kripke/profiler.h"
#endif
#include "Kripke/transport_headers.h"

int main(int argc, char *argv[])
{
  FILE             *in_file;
  FILE             *out_file;
  Input_Variables  *input_variables;
  int myid;
  int ierr=0;
  int timing_index[3];
  int i;
  char input_file_name[256];
  char output_file_name[256];

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
  else {
    strcpy(input_file_name, argv[1]);
  }
  if((in_file = fopen(input_file_name, "r")) == NULL){
    printf("ERROR: Can't open input file %s.\n", input_file_name);
    error_exit(1);
  }
  if(myid == 0){
    /* Print out message about input file name. */
    printf("\n");
    printf("Input File: %s\n", input_file_name);
    printf("\n");
  }

  /*-----------------------------------------------------------------------
   * Initialize timing
   *-----------------------------------------------------------------------*/
  timing_index[0] = InitializeTiming("Total Run Time");
  timing_index[1] = InitializeTiming("Setup");
  timing_index[2] = InitializeTiming("Total Solve");
  BeginTiming(timing_index[0]);

  /*-----------------------------------------------------------------------
   * Initialize user input and put relevant data in the user_data structure
   *-----------------------------------------------------------------------*/
  BeginTiming(timing_index[1]);
  input_variables = ReadInput(in_file);
  if(myid == 0){
    sprintf(output_file_name, "%s.out", input_variables->run_name);
    if((out_file = fopen(output_file_name, "w")) == NULL){
      printf("ERROR: Can't open output file %s.\n", output_file_name);
      error_exit(1);
    }
    PrintInputVariables(input_variables, out_file);
    fclose(out_file);
  }
  User_Data *user_data = new User_Data(MPI_COMM_WORLD, input_variables);
  FreeInputVariables(input_variables);


  EndTiming(timing_index[1]);

  /* Begin timing of solve */
  BeginTiming(timing_index[2]);

  /* Run the driver */
  for(i=0; i<user_data->ncalls; i++){
    SweepDriver(user_data);
  }

  /* End timing of solve */
  EndTiming(timing_index[2]);

  /* End timing of total run time */
  EndTiming(timing_index[0]);

  /*-----------------------------------------------------------------------
   * Print timing information
   *-----------------------------------------------------------------------*/
  PrintTiming("SWEEP KERNEL Run Time Statistics", MPI_COMM_WORLD);

  /* Free timing information */
  FinalizeTiming(timing_index[0]);
  FinalizeTiming(timing_index[1]);
  FinalizeTiming(timing_index[2]);
  for(i=0; i<user_data->num_timings; i++){
    FinalizeTiming(user_data->timing_index[i]);
  }
  FREE(user_data->timing_index);
  delete user_data;

  MPI_Finalize();

  return(0);
}
