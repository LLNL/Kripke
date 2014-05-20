/*--------------------------------------------------------------------------
 * Main program for 3-D Sweep Kernel
 *--------------------------------------------------------------------------*/
#include <Kripke.h>
#include <Kripke/User_Data.h>
#include <Kripke/Comm.h>

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<sstream>

#ifdef KRIPKE_USE_OPENMP
//#include<openmp.h>
#endif

#ifdef KRIPKE_USE_PERFTOOLS
#include<google/profiler.h>
#endif

int main(int argc, char *argv[])
{
  FILE             *in_file;
  FILE             *out_file;
  int myid;
  int ierr=0;
  int i;

  bool profile = true;

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
    printf("-------------- KRIPKE SWEEPER VERSION 1.0 ---------------\n");
    printf("---------------------------------------------------------\n");
#ifdef KRIPKE_USE_OPENMP
  printf("Using OpenMP\n");
#endif
  }

  if(argc < 2){
    printf("ERROR: Please specify an input file name.\n");
    error_exit(1);
  }
  std::string input_file_name = argv[1];


  /*-----------------------------------------------------------------------
   * Initialize user input and put relevant data in the user_data structure
   *-----------------------------------------------------------------------*/
#ifdef KRIPKE_USE_PERFTOOLS
  if(profile){
    std::stringstream pfname;
    pfname << "profile." << myid;
    ProfilerStart(pfname.str().c_str());
    ProfilerRegisterThread();
  }
#endif



  Input_Variables input_variables;
  input_variables.read(input_file_name);
  input_variables.print();
  User_Data *user_data = new User_Data(&input_variables);

  std::vector<std::string> papi_names;
  for(int i = 2;i < argc; ++i){
    papi_names.push_back(argv[i]);
    //papi_names.push_back("PAPI_TOT_CYC");
    //papi_names.push_back("PAPI_FP_INS");
    //papi_names.push_back("PAPI_VEC_DP");
  }
  user_data->timing.setPapiEvents(papi_names);



  /* Run the driver */
  Driver(user_data);

#ifdef KRIPKE_USE_PERFTOOLS
  if(profile){
    ProfilerFlush();
    ProfilerStop();
  }
#endif


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
