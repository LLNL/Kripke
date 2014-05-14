#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#include <Kripke.h>
#include <Kripke/User_Data.h>
#include <Kripke/Comm.h>

struct runLTimes {
  std::string name(void) const { return "LTimes"; }

  void operator ()(User_Data *user_data) const {
    user_data->kernel->LTimes(user_data->grid_data);
  }
};

struct runLPlusTimes {
  std::string name(void) const { return "LPlusTimes"; }

  void operator ()(User_Data *user_data) const {
    user_data->kernel->LPlusTimes(user_data->grid_data);
  }
};

struct runSweep {
  std::string name(void) const { return "Sweep"; }

  void operator ()(User_Data *user_data) const {
    for(int group_set = 0;group_set < user_data->num_group_sets;++ group_set){
      SweepSolver_GroupSet(group_set, user_data);
    }
  }
};


template<typename KernelRunner>
void testKernel(Input_Variables &input_variables, Nesting_Order nesting){
  KernelRunner kr;

  printf("Comparing %s to %s for kernel %s\n",
      nestingString(NEST_GDZ).c_str(),
      nestingString(nesting).c_str(),
      kr.name().c_str());

  // Allocate two problems (one reference)
  printf(" -- allocating\n");
  input_variables.nesting = nesting;
  User_Data *user_data = new User_Data(&input_variables);

  input_variables.nesting = NEST_GDZ;
  User_Data *ref_data = new User_Data(&input_variables);

  // Generate random data in the reference problem, and copy it to the other
  printf(" -- randomizing data\n");
  ref_data->randomizeData();
  user_data->copy(*ref_data);

  printf(" -- running kernels\n");

  // Run both kernels
  kr(ref_data);
  kr(user_data);

  printf(" -- comparing results\n");
  // Compare differences
  bool is_diff = ref_data->compare(*user_data, 1e-12, true);
  if(is_diff){
    printf("Differences found, bailing out\n");
    error_exit(1);
  }

  // Cleanup
  printf(" -- done\n");
  delete user_data;
  delete ref_data;
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
    printf("------------------- KRIPKE KERNEL TEST ------------------\n");
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
/*
  // Run LTimes
  for(int nesting = 0;nesting < 6; ++ nesting){
    testKernel<runLTimes>(input_variables, (Nesting_Order)nesting);
  }

  // Run LPlusTimes
  for(int nesting = 0;nesting < 6; ++ nesting){
    testKernel<runLPlusTimes>(input_variables, (Nesting_Order)nesting);
  } */

  // Run Sweep
  for(int nesting = 0;nesting < 6; ++ nesting){
    testKernel<runSweep>(input_variables, (Nesting_Order)nesting);
  }

  /*-----------------------------------------------------------------------
   * Cleanup and exist
   *-----------------------------------------------------------------------*/
  MPI_Finalize();
  return(0);
}
