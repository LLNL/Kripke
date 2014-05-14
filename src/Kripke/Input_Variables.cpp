/*--------------------------------------------------------------------------
 * Input-reading utilities
 *--------------------------------------------------------------------------*/

#include <Kripke/Input_Variables.h>
#include <Kripke/Comm.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <map>
#include <mpi.h>

void Input_Variables::read(std::string const &fname)
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  // if we are the root processor, read in the file
  if(rank == 0){
    FILE *in_file;
    printf("\nInput File: %s\n\n", fname.c_str());
    in_file = fopen(fname.c_str(), "rb");
    if(in_file == NULL){
      printf("Couldn't open input file\n");
      error_exit(1);
    }

    /*
     * Define mappings from variable names to our struct
     */
    std::map<std::string, int *> int_vars;

    int_vars["npx"] = &npx;
    int_vars["npy"] = &npy;
    int_vars["npz"] = &npz;
    int_vars["niter"] = &niter;
    int_vars["nx"] = &nx;
    int_vars["ny"] = &ny;
    int_vars["nz"] = &nz;

    int_vars["num_dirsets_per_octant"] = &num_dirsets_per_octant;
    int_vars["num_dirs_per_dirset"] = &num_dirs_per_dirset;
    int_vars["num_groupsets"] = &num_groupsets;
    int_vars["num_groups_per_groupset"] = &num_groups_per_groupset;

    int_vars["nesting"] = (int*)&nesting;
    int_vars["sweep_order"] = (int*)&sweep_order;
    int_vars["block_size"] = &block_size;

    /*
     * Parse input file
     */
    while(!feof(in_file) ){
      char line[1024];
      fgets(line, 1024, in_file);

      // strip white space
      char *line_start;
      for(line_start = line; *line_start != 0; ++line_start){
        if(!isspace(*line_start)){
          break;
        }
      }

      // strip comments
      for(char *ptr = line_start; *ptr != 0;++ptr){
        if(*ptr == '#'){
          *ptr = 0;
          break;
        }
      }

      // anything left?
      if(strlen(line_start) == 0){
        continue;
      }

      // parse it
      char key[1024];
      int val;
      sscanf(line_start, "%s%d", key, &val);
      if(int_vars.find(key) != int_vars.end()){
        (*int_vars[key]) = val;
      }
      else{
        printf("Unknown input: '%s'\n", line);
        error_exit(1);
      }
    }
  }

  // broadcast the data to everyone
  MPI_Bcast(this, sizeof(Input_Variables), MPI_BYTE, 0, MPI_COMM_WORLD);
}

/*--------------------------------------------------------------------------
 * PrintInputVariables : Prints values in an Input_Variables structure
 *--------------------------------------------------------------------------*/

void Input_Variables::print(void) const
/*--------------------------------------------------------------------------
 * input_variables : The Input_Variables structure to be printed.
 * out_file        : The file to which input information will be printed.
 *--------------------------------------------------------------------------*/
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if(rank != 0){
    return;
  }
  printf("\n");
  printf("Processors\n");
  printf("   npx = %i   npy = %i   npz = %i\n",
          npx, npy, npz);
  printf("   niter = %i\n", niter);
  printf("\n");
  printf("Space_Grid\n");
  printf("   nx = %i  ny = %i  nz = %i\n",
          nx, ny, nz);
  printf("\n");
  printf("DiscreteOrdinates_Grid\n");
  printf("   num_dirsets_per_octant =   %d\n", num_dirsets_per_octant);
  printf("   num_dirs_per_dirset =      %d\n", num_dirs_per_dirset);
  printf("   total directions =         %d\n", 8*num_dirsets_per_octant*num_dirs_per_dirset);
  printf("\n");
  printf("Energy_Grid\n");
  printf("   num_groupsets =           %d\n", num_groupsets);
  printf("   num_groups_per_groupset = %d\n", num_groups_per_groupset);
  printf("   total groups =            %d\n", num_groupsets * num_groups_per_groupset);
  printf("\n");
  printf("Nesting: GS - DS - ");
  switch(nesting){
    case NEST_DGZ:
      printf("D - G - Z\n");
      break;
    case NEST_GDZ:
      printf("G - D - Z\n");
      break;
    case NEST_DZG:
      printf("D - Z - G\n");
      break;
    case NEST_ZDG:
      printf("Z - D - G\n");
      break;
    case NEST_GZD:
      printf("G - Z - D\n");
      break;
    case NEST_ZGD:
      printf("Z - G - D\n");
      break;
    default:
      printf(" UNKNOWN\n");
      error_exit(1);
      break;
  }

  printf("\n");
}

