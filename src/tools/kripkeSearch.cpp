#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#include <Kripke.h>
#include <Kripke/User_Data.h>
#include <Kripke/Comm.h>

void runPoint(Input_Variables &input_variables, FILE *out_fp){
  /* Allocate problem */
  User_Data *user_data = new User_Data(&input_variables);

  /* Run the driver */
  Driver(user_data);

  /* Get timing data */
  //user_data->timing.print();

  std::string nesting = nestingString(input_variables.nesting);

  char line[2048];
  double niter = (double)input_variables.niter;
  snprintf(line, 2048, "%d %s %3d %3d %3d %3d %8.4lf %8.4lf %8.4lf %8.4lf\n",
      (int)input_variables.nesting,
      nesting.c_str(),
      input_variables.num_dirsets_per_octant,
      input_variables.num_dirs_per_dirset,
      input_variables.num_groupsets,
      input_variables.num_groups_per_groupset,
      user_data->timing.getTotal("Solve")/niter,
      user_data->timing.getTotal("Sweep")/niter,
      user_data->timing.getTotal("LTimes")/niter,
      user_data->timing.getTotal("LPlusTimes")/niter
    );
  if(out_fp != NULL){
    fprintf(out_fp, line);
    fflush(out_fp);
    printf(line);
  }

  /* Cleanup */
  delete user_data;
}

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
    printf("--------------- KRIPKE SEARCH VERSION 1.0 ---------------\n");
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


  FILE *fp = NULL;
  if(myid == 0){
    std::string output_file_name = input_file_name + ".out";
    fp = fopen(output_file_name.c_str(), "wb");
  }

  int total_dpow = 4;
  int total_gpow = 4;

  for(int dset_pow = 0;dset_pow <= total_dpow;++dset_pow){
    int n_dset = 1 << dset_pow;
    int n_dirs = 1 << (total_dpow-dset_pow);

    for(int gset_pow = 0;gset_pow <= total_gpow;++gset_pow){
      int n_gset = 1 << gset_pow;
      int n_grps = 1 << (total_gpow-gset_pow);

      for(int nesting = 0;nesting < 6; ++ nesting){
        input_variables.nesting = (Nesting_Order)nesting;
        input_variables.num_dirsets_per_octant = n_dset;
        input_variables.num_dirs_per_dirset = n_dirs;
        input_variables.num_groupsets = n_gset;
        input_variables.num_groups_per_groupset = n_grps;
        runPoint(input_variables, fp);
      }
    }
  }

  if(fp != NULL){
    fclose(fp);
  }

  /*-----------------------------------------------------------------------
   * Cleanup and exist
   *-----------------------------------------------------------------------*/
  MPI_Finalize();
  return(0);
}
