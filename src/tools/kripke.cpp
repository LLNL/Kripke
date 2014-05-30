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
    printf("  --procs <npx,npy,npz>  MPI task spatial decomposition\n");
    printf("                         Default:  --procs 1,1,1\n");
    printf("  --test                 Run Kernel Test instead of solver\n");
    printf("  --zones <x,y,z>        Number of zones in x,y,z\n");
    printf("                         Default:  --zones 8,8,8\n");
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


void runPoint(Input_Variables &input_variables, FILE *out_fp){
  /* Allocate problem */
  User_Data *user_data = new User_Data(&input_variables);

  user_data->timing.setPapiEvents(papi_names);

  /* Run the driver */
  Driver(user_data);

  std::string nesting = nestingString(input_variables.nesting);

  char line[2048];
  double niter = (double)input_variables.niter;
  snprintf(line, 2048, "nestid=%d nest=%s D=%-3d d=%-3d dirs=%d G=%-3d g=%-3d grps=%d Solve=%-8.4lf Sweep=%-8.4lf LTimes=%-8.4lf LPlusTimes=%-8.4lf\n",
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
  if(out_fp != NULL){
    user_data->timing.print();

    if(out_fp != stdout){
      fprintf(out_fp, line);
    }
    printf(line);
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
  int nzones[3] = {8, 8, 8};
  int lorder = 2;
  int niter = 10;
  bool test = false;

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
    else{
      printf("Unknwon options %s\n", opt.c_str());
      usage();
    }
  }

  /*
   * Display Options
   */
  int nsearches = grp_list.size() * dir_list.size() * nest_list.size();
  if (myid == 0) {
    printf("Number of MPI tasks:   %d\n", num_tasks);
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
  }

  /*
   * Execute the Search Space
   */
  FILE *outfp = NULL;
  if(outfile != "" && myid == 0){
    outfp = fopen(outfile.c_str(), "wb");
  }
  else if(myid == 0){
    outfp = stdout;
  }
  Input_Variables ivars;
  ivars.nx = nzones[0];
  ivars.ny = nzones[1];
  ivars.nz = nzones[2];
  ivars.npx = nprocs[0];
  ivars.npy = nprocs[1];
  ivars.npz = nprocs[2];
  ivars.block_size = 0;
  ivars.legendre_order = lorder + 1;
  ivars.sweep_order = SWEEP_DEFAULT;
  ivars.niter = niter;
  int point = 0;
  for(int d = 0;d < dir_list.size();++ d){
    for(int g = 0;g < grp_list.size();++ g){
      for(int n = 0;n < nest_list.size();++ n){

        printf("Running point %d/%d: D:d=%d:%d, G:g=%d:%d, Nest=%s\n",
            point+1, nsearches,
            dir_list[d].first,
            dir_list[d].second,
            grp_list[d].first,
            grp_list[d].second,
            nestingString(nest_list[n]).c_str());

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
          runPoint(ivars, outfp);
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
  return (0);
}
