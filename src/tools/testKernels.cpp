#include "testKernels.h"

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
void testKernel(Input_Variables &input_variables){
  KernelRunner kr;

  printf("  Comparing %s to %s for kernel %s\n",
      nestingString(NEST_GDZ).c_str(),
      nestingString(input_variables.nesting).c_str(),
      kr.name().c_str());

  // Allocate two problems (one reference)
  printf("    -- allocating\n");
  User_Data *user_data = new User_Data(&input_variables);

  Nesting_Order old_nest = input_variables.nesting;
  input_variables.nesting = NEST_GDZ;
  User_Data *ref_data = new User_Data(&input_variables);
  input_variables.nesting = old_nest;

  // Generate random data in the reference problem, and copy it to the other
  printf("    -- randomizing data\n");
  ref_data->randomizeData();
  user_data->copy(*ref_data);

  printf("    -- running kernels\n");

  // Run both kernels
  kr(ref_data);
  kr(user_data);

  printf("    -- comparing results\n");
  // Compare differences
  bool is_diff = ref_data->compare(*user_data, 1e-12, true);
  if(is_diff){
    printf("Differences found, bailing out\n");
    error_exit(1);
  }

  // Cleanup
  printf("    -- OK\n\n");
  delete user_data;
  delete ref_data;
}


void testKernels(Input_Variables &input_variables){
  // Run LTimes
  testKernel<runLTimes>(input_variables);

  // Run LPlusTimes
  testKernel<runLPlusTimes>(input_variables);

  // Run Sweep
  testKernel<runSweep>(input_variables);
}
