/**
 * This file contains all of the correctness checking code for Kripke.
 */
#include "testKernels.h"

#include <Kripke.h>
#include <Kripke/User_Data.h>
#include <Kripke/Comm.h>

/**
 * Functional object to run the LTimes kernel.
 */
struct runLTimes {
  std::string name(void) const { return "LTimes"; }

  void operator ()(User_Data *user_data) const {
    user_data->kernel->LTimes(user_data->grid_data);
  }
};

/**
 * Functional object to run the LPlusTimes kernel.
 */
struct runLPlusTimes {
  std::string name(void) const { return "LPlusTimes"; }

  void operator ()(User_Data *user_data) const {
    user_data->kernel->LPlusTimes(user_data->grid_data);
  }
};

/**
 * Functional object to run the MPI sweep and sweep kernels
 */
struct runSweep {
  std::string name(void) const { return "Sweep"; }

  void operator ()(User_Data *user_data) const {
    for(int group_set = 0;group_set < user_data->num_group_sets;++ group_set){
      SweepSolver_GroupSet(group_set, user_data);
    }
  }
};


/**
 * Tests a specific kernel (using one of the above runXXX functional objects).
 */
template<typename KernelRunner>
void testKernel(Input_Variables &input_variables){
  int myid;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  KernelRunner kr;

  if(myid == 0){
    printf("  Comparing %s to %s for kernel %s\n",
      nestingString(NEST_GDZ).c_str(),
      nestingString(input_variables.nesting).c_str(),
      kr.name().c_str());
  }

  // Allocate two problems (one reference)
  if(myid == 0)printf("    -- allocating\n");
  User_Data *user_data = new User_Data(&input_variables);

  Nesting_Order old_nest = input_variables.nesting;
  input_variables.nesting = NEST_GDZ;
  User_Data *ref_data = new User_Data(&input_variables);
  input_variables.nesting = old_nest;

  // Generate random data in the reference problem, and copy it to the other
  if(myid == 0)printf("    -- randomizing data\n");
  ref_data->randomizeData();
  user_data->copy(*ref_data);

  if(myid == 0)printf("    -- running kernels\n");

  // Run both kernels
  kr(ref_data);
  kr(user_data);

  if(myid == 0)printf("    -- comparing results\n");
  // Compare differences
  bool is_diff = ref_data->compare(*user_data, 1e-12, true);
  if(is_diff){
    if(myid == 0)printf("Differences found, bailing out\n");
    error_exit(1);
  }

  // Cleanup
  if(myid == 0)printf("    -- OK\n\n");
  delete user_data;
  delete ref_data;
}


/**
 * Tests all kernels given the specified input.
 */
void testKernels(Input_Variables &input_variables){
  // Run LTimes
  testKernel<runLTimes>(input_variables);

  // Run LPlusTimes
  testKernel<runLPlusTimes>(input_variables);

  // Run Sweep
  testKernel<runSweep>(input_variables);
}
