#include<Kripke.h>
#include<Kripke/Input_Variables.h>
#include<Kripke/User_Data.h>
#include"testKernels.h"
#include<stdio.h>
#include<string.h>
#include<mpi.h>
#include<algorithm>
#include<string>
#include<sstream>

#ifdef KRIPKE_USE_OPENMP
#include<omp.h>
#endif

#ifdef KRIPKE_USE_TCMALLOC
#include<gperftools/malloc_extension.h>
#endif

#ifdef KRIPKE_USE_PERFTOOLS
#include<google/profiler.h>
#endif

#ifdef __bgq__
#include </bgsys/drivers/ppcfloor/spi/include/kernel/location.h>
#include </bgsys/drivers/ppcfloor/spi/include/kernel/memory.h>
#endif

typedef std::pair<int, int> IntPair;

std::vector<std::string> papi_names;

void usage(void){
  int myid;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  if(myid == 0){
    printf("Usage:  [srun ...] kripke [options...]\n");
    printf("Where options are:\n");
    printf("  --dir [D:d,D:d,...]    List of dirsets and dirs/set pairs\n");
    printf("                         Default:  --dir 1:1\n");
    printf("                         Example:  --dir 1:4,2:2,4:1\n");
    printf("  --grp [G:g,G:g,...]    List of grpsets and groups/set pairs\n");
    printf("                         Default:  --grp 1:1\n");
    printf("  --legendre <lorder>    Scattering Legendre Expansion Order (0, 1, ...)\n");
    printf("                         Default:  --legendre 2\n");
    printf("  --nest [n,n,...]       List of data nestings\n");
    printf("                         Default:  --nest DGZ,DZG,GDZ,GZD,ZDG,ZGD\n");
    printf("  --niter <NITER>        Number of solver iterations to run (default: 10)\n");
    printf("  --out <OUTFILE>        Optional output file (default: none)\n");
    printf("  --gperf                Turn on Google Perftools profiling\n");
    printf("  --procs <npx,npy,npz>  MPI task spatial decomposition\n");
    printf("                         Default:  --procs 1,1,1\n");
    printf("  --restart <point>      Restart at given point\n");
    printf("  --test                 Run Kernel Test instead of solver\n");
    printf("  --zones <x,y,z>        Number of zones in x,y,z\n");
    printf("                         Default:  --zones 12,12,12\n");
    printf("\n");
  }
  MPI_Finalize();
  exit(1);
}

struct CmdLine {
  CmdLine(int argc, char **argv) :
    size(argc-1),
    cur(0),
    args()
  {
    for(int i = 0;i < size;++ i){
      args.push_back(argv[i+1]);
    }
  }

  std::string pop(void){
    if(atEnd())
      usage();
    return args[cur++];
  }

  bool atEnd(void){
    return(cur >= size);
  }

  int size;
  int cur;
  std::vector<std::string> args;
};

std::vector<std::string> split(std::string const &str, char delim){
  std::vector<std::string> elem;
  std::stringstream ss(str);
  std::string e;
  while(std::getline(ss, e, delim)){
    elem.push_back(e);
  }
  return elem;
}


void runPoint(int point, int num_tasks, int num_threads, Input_Variables &input_variables, FILE *out_fp){

  /* Allocate problem */
  User_Data *user_data = new User_Data(&input_variables);

  user_data->timing.setPapiEvents(papi_names);

  /* Run the solver */
  SweepSolver(user_data);

  std::string nesting = nestingString(input_variables.nesting);

  char line[2048];
  double niter = (double)input_variables.niter;
  snprintf(line, 2048, "RUN: point=%d ntasks=%d nthreads=%d nestid=%d nest=%s D=%-3d d=%-3d dirs=%d G=%-3d g=%-3d grps=%d Solve=%-8.4lf Sweep=%-8.4lf LTimes=%-8.4lf LPlusTimes=%-8.4lf\n",
      point,
      num_tasks,
      num_threads,
      (int)input_variables.nesting,
      nesting.c_str(),
      input_variables.num_dirsets_per_octant,
      input_variables.num_dirs_per_dirset,
      8*input_variables.num_dirsets_per_octant*input_variables.num_dirs_per_dirset,
      input_variables.num_groupsets,
      input_variables.num_groups_per_groupset,
      input_variables.num_groupsets * input_variables.num_groups_per_groupset,
      user_data->timing.getTotal("Solve")/niter,
      user_data->timing.getTotal("Sweep")/niter,
      user_data->timing.getTotal("LTimes")/niter,
      user_data->timing.getTotal("LPlusTimes")/niter
    );
  int myid;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  if(myid == 0){
    if(out_fp != NULL){
      fprintf(out_fp, line);
      fflush(out_fp);
    }
    user_data->timing.print();
    printf(line);
    printf("\n\n");
  }

  /* Cleanup */
  delete user_data;
}

int main(int argc, char **argv) {
  /*
   * Initialize MPI
   */
  MPI_Init(&argc, &argv);
  int myid;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  int num_tasks;
  MPI_Comm_size(MPI_COMM_WORLD, &num_tasks);

  if (myid == 0) {
    /* Print out a banner message along with a version number. */
    printf("\n");
    printf("---------------------------------------------------------\n");
    printf("------------------- KRIPKE VERSION 1.0 ------------------\n");
    printf("---------------------------------------------------------\n");

    /* Print out some information about how OpenMP threads are being mapped
     * to CPU cores.
     */
#ifdef KRIPKE_USE_OPENMP
#pragma omp parallel
    {
      int tid = omp_get_thread_num();
#ifdef __bgq__
      int core = Kernel_ProcessorCoreID();
#else
      int core = sched_getcpu();
#endif
      printf("Rank: %d Thread %d: Core %d\n", myid, tid, core);
    }
#endif
  }


  /*
   * Default input parameters
   */
  std::vector<IntPair> grp_list;
  grp_list.push_back(IntPair(1,1));
  std::vector<IntPair> dir_list;
  dir_list.push_back(IntPair(1,1));
  std::string outfile;
  int nprocs[3] = {1, 1, 1};
  int nzones[3] = {12, 12, 12};
  int lorder = 4;
  int niter = 10;
  bool test = false;
  bool perf_tools = false;
  int restart_point = 0;

  std::vector<Nesting_Order> nest_list;
  nest_list.push_back(NEST_DGZ);
  nest_list.push_back(NEST_DZG);
  nest_list.push_back(NEST_GDZ);
  nest_list.push_back(NEST_GZD);
  nest_list.push_back(NEST_ZDG);
  nest_list.push_back(NEST_ZGD);

  /*
   * Parse command line
   */
  CmdLine cmd(argc, argv);
  while(!cmd.atEnd()){
    std::string opt = cmd.pop();
    if(opt == "-h" || opt == "--help"){usage();}
    else if(opt == "--out"){outfile = cmd.pop();}
    else if(opt == "--zones"){
      std::vector<std::string> nz = split(cmd.pop(), ',');
      if(nz.size() != 3) usage();
      nzones[0] = std::atoi(nz[0].c_str());
      nzones[1] = std::atoi(nz[1].c_str());
      nzones[2] = std::atoi(nz[2].c_str());
    }
    else if(opt == "--procs"){
      std::vector<std::string> np = split(cmd.pop(), ',');
      if(np.size() != 3) usage();
      nprocs[0] = std::atoi(np[0].c_str());
      nprocs[1] = std::atoi(np[1].c_str());
      nprocs[2] = std::atoi(np[2].c_str());
    }
    else if(opt == "--grp"){
      std::vector<std::string> sets = split(cmd.pop(), ',');
      if(sets.size() < 1)usage();
      grp_list.clear();
      for(int i = 0;i < sets.size();++ i){
        std::vector<std::string> p = split(sets[i], ':');
        if(p.size() != 2)usage();
        grp_list.push_back(IntPair(std::atoi(p[0].c_str()), std::atoi(p[1].c_str())));
      }
    }
    else if(opt == "--dir"){
      std::vector<std::string> sets = split(cmd.pop(), ',');
      if(sets.size() < 1)usage();
      dir_list.clear();
      for(int i = 0;i < sets.size();++ i){
        std::vector<std::string> p = split(sets[i], ':');
        if(p.size() != 2)usage();
        dir_list.push_back(IntPair(std::atoi(p[0].c_str()), std::atoi(p[1].c_str())));
      }
    }
    else if(opt == "--legendre"){
      lorder = std::atoi(cmd.pop().c_str());
    }
    else if(opt == "--niter"){
      niter = std::atoi(cmd.pop().c_str());
    }
    else if(opt == "--nest"){
      std::vector<std::string> sets = split(cmd.pop(), ',');
      if(sets.size() < 1)usage();
      nest_list.clear();
      for(int i = 0;i < sets.size();++ i){
        Nesting_Order n = nestingFromString(sets[i]);
        if(n < 0)usage();
        nest_list.push_back(n);
      }
    }
    else if(opt == "--test"){
      test = true;
    }
    else if(opt == "--papi"){
      papi_names = split(cmd.pop(), ',');
    }
    else if(opt == "--gperf"){
      perf_tools = true;
    }
    else if(opt == "--restart"){
      restart_point = std::atoi(cmd.pop().c_str());
    }
    else{
      printf("Unknwon options %s\n", opt.c_str());
      usage();
    }
  }

  /*
   * Display Options
   */
  int nsearches = grp_list.size() * dir_list.size() * nest_list.size();
  int num_threads=1;
  if (myid == 0) {
    printf("Number of MPI tasks:   %d\n", num_tasks);
#ifdef KRIPKE_USE_OPENMP
#pragma omp parallel
    {
      num_threads = omp_get_num_threads();
      if(omp_get_thread_num() == 0){
          printf("OpenMP threads/task:   %d\n", num_threads);
          printf("OpenMP total threads:  %d\n", num_threads*num_tasks);
        }
    }
#endif
    printf("Output File:           %s\n", outfile.c_str());
    printf("Processors:            %d x %d x %d\n", nprocs[0], nprocs[1], nprocs[2]);
    printf("Zones:                 %d x %d x %d\n", nzones[0], nzones[1], nzones[2]);
    printf("Legendre Order:        %d\n", lorder);
    printf("Number iterations:     %d\n", niter);

    if(grp_list.size() == 0){
      printf("No GroupSet/Groups defined (--grp)\n");
      usage();
    }
    printf("GroupSet/Groups:       ");
    for(int i = 0;i < grp_list.size();++ i){
      printf("%s%d:%d", (i==0 ? "" : ", "), grp_list[i].first, grp_list[i].second);
    }
    printf("\n");

    if(dir_list.size() == 0){
      printf("No DirSets/Directions defined (--dir)\n");
      usage();
    }
    printf("DirSets/Directions:    ");
    for(int i = 0;i < dir_list.size();++ i){
      printf("%s%d:%d", (i==0 ? "" : ", "), dir_list[i].first, dir_list[i].second);
    }
    printf("\n");

    printf("Nestings:              ");
    for(int i = 0;i < nest_list.size();++ i){
      printf("%s%s", (i==0 ? "" : ", "), nestingString(nest_list[i]).c_str());
    }
    printf("\n");
    printf("Search space size:     %d points\n", nsearches);
    if(perf_tools){
      printf("Using Google Perftools\n");
    }
  }

  /*
   * Execute the Search Space
   */
  FILE *outfp = NULL;
  if(outfile != "" && myid == 0){
    if(restart_point == 0){
      outfp = fopen(outfile.c_str(), "wb");
    }
    else{
      outfp = fopen(outfile.c_str(), "ab");
    }
  }
#ifdef KRIPKE_USE_PERFTOOLS
  if(perf_tools){
    std::stringstream pfname;
    pfname << "profile." << myid;
    ProfilerStart(pfname.str().c_str());
    ProfilerRegisterThread();
  }
#endif
  Input_Variables ivars;
  ivars.nx = nzones[0];
  ivars.ny = nzones[1];
  ivars.nz = nzones[2];
  ivars.npx = nprocs[0];
  ivars.npy = nprocs[1];
  ivars.npz = nprocs[2];
  ivars.legendre_order = lorder + 1;
  ivars.block_size = 4;
  ivars.niter = niter;
  int point = 0;
  for(int d = 0;d < dir_list.size();++ d){
    for(int g = 0;g < grp_list.size();++ g){
      for(int n = 0;n < nest_list.size();++ n){

        if(restart_point <= point+1){
          if(myid == 0){
            printf("Running point %d/%d: D:d=%d:%d, G:g=%d:%d, Nest=%s\n",
                point+1, nsearches,
                dir_list[d].first,
                dir_list[d].second,
                grp_list[g].first,
                grp_list[g].second,
                nestingString(nest_list[n]).c_str());
          }
          // Setup Current Search Point
          ivars.num_dirsets_per_octant = dir_list[d].first;
          ivars.num_dirs_per_dirset = dir_list[d].second;
          ivars.num_groupsets = grp_list[g].first;
          ivars.num_groups_per_groupset = grp_list[g].second;
          ivars.nesting = nest_list[n];

          // Run the point
          if(test){
            // Invoke Kernel testing
            testKernels(ivars);
          }
          else{
            // Just run the "solver"
            runPoint(point+1, num_tasks, num_threads, ivars, outfp);
          }


          // Gather post-point memory info
          double heap_mb = -1.0;
          double hwm_mb = -1.0;
#ifdef KRIPKE_USE_TCMALLOC
          // If we are using tcmalloc, we need to use it's interface
          MallocExtension *mext = MallocExtension::instance();
          size_t bytes;

          mext->GetNumericProperty("generic.current_allocated_bytes", &bytes);
          heap_mb = ((double)bytes)/1024.0/1024.0;

          mext->GetNumericProperty("generic.heap_size", &bytes);
          hwm_mb = ((double)bytes)/1024.0/1024.0;
#else
#ifdef __bgq__
          // use BG/Q specific calls (if NOT using tcmalloc)
          uint64_t bytes;

          rc = Kernel_GetMemorySize(KERNEL_MEMSIZE_HEAP, &bytes);
          heap_mb = ((double)bytes)/1024.0/1024.0;

          rc = Kernel_GetMemorySize(KERNEL_MEMSIZE_HEAPMAX, &bytes);
          hwm_mb = ((double)bytes)/1024.0/1024.0;
#endif
#endif
          // Print memory info
          if(myid == 0 && heap_mb >= 0.0){
            printf("Bytes allocated: %lf MB\n", heap_mb);
            printf("Heap Size      : %lf MB\n", hwm_mb);

          }
        }
        point ++;

      }
    }
  }
  if(outfp != NULL){
    fclose(outfp);
  }

  /*
   * Cleanup and exit
   */
  MPI_Finalize();
#ifdef KRIPKE_USE_PERFTOOLS
  if(perf_tools){
    ProfilerFlush();
    ProfilerStop();
  }
#endif
  return (0);
}
