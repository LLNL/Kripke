/*--------------------------------------------------------------------------
 * Input-reading utilities
 *--------------------------------------------------------------------------*/

#include <Kripke/input_variables.h>
#include <Kripke/comm.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>



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

    char file_line[256];
    char key_name[256];
    char test_key[256];
    char str_tmp[256];

    int n, i, j;

    int int_tmp;
    double double_tmp;

    while(!feof(in_file) ){
      fgets(file_line, 256, in_file);

      if(file_line[0] != '#'){
        sscanf(file_line, "%s", key_name);

        if(strcmp(key_name, "Run_Name") == 0){
          fgets(file_line, 256, in_file);
          while(file_line[0] == '#'){
            fgets(file_line, 256, in_file);
          }
          sscanf(file_line, "%s %s", key_name, str_tmp);
          if(strcmp(key_name, "run_name") == 0){
            strcpy(run_name, str_tmp);
          }
        }
        else if(strcmp(key_name, "Processors") == 0){
          fgets(file_line, 256, in_file);
          while(file_line[0] == '#'){
            fgets(file_line, 256, in_file);
          }
          sscanf(file_line, "%s %d", key_name, &int_tmp);
          if(strcmp(key_name, "npx") == 0){
            npx = int_tmp;
          }
          else {printf("Error reading npx.\n"); }

          fgets(file_line, 256, in_file);
          while(file_line[0] == '#'){
            fgets(file_line, 256, in_file);
          }
          sscanf(file_line, "%s %d", key_name, &int_tmp);
          if(strcmp(key_name, "npy") == 0){
            npy = int_tmp;
          }
          else {printf("Error reading npy.\n"); }

          fgets(file_line, 256, in_file);
          while(file_line[0] == '#'){
            fgets(file_line, 256, in_file);
          }
          sscanf(file_line, "%s %d", key_name, &int_tmp);
          if(strcmp(key_name, "npz") == 0){
            npz = int_tmp;
          }
          else {printf("Error reading npz.\n"); }

          fgets(file_line, 256, in_file);
          while(file_line[0] == '#'){
            fgets(file_line, 256, in_file);
          }
          sscanf(file_line, "%s %d", key_name, &int_tmp);
          if(strcmp(key_name, "ncalls") == 0){
            ncalls = int_tmp;
          }
          else {printf("Error reading nlevels_kba.\n"); }
        }
        else if(strcmp(key_name, "Domain") == 0){
          fgets(file_line, 256, in_file);
          while(file_line[0] == '#'){
            fgets(file_line, 256, in_file);
          }
          sscanf(file_line, "%s %lf", key_name, &double_tmp);
          if(strcmp(key_name, "xmin") == 0){
            xmin = double_tmp;
          }
          else {printf("Error reading xmin.\n"); }

          fgets(file_line, 256, in_file);
          while(file_line[0] == '#'){
            fgets(file_line, 256, in_file);
          }
          sscanf(file_line, "%s %lf", key_name, &double_tmp);
          if(strcmp(key_name, "xmax") == 0){
            xmax = double_tmp;
          }
          else {printf("Error reading xmax.\n"); }

          fgets(file_line, 256, in_file);
          while(file_line[0] == '#'){
            fgets(file_line, 256, in_file);
          }
          sscanf(file_line, "%s %lf", key_name, &double_tmp);
          if(strcmp(key_name, "ymin") == 0){
            ymin = double_tmp;
          }
          else {printf("Error reading ymin.\n"); }

          fgets(file_line, 256, in_file);
          while(file_line[0] == '#'){
            fgets(file_line, 256, in_file);
          }
          sscanf(file_line, "%s %lf", key_name, &double_tmp);
          if(strcmp(key_name, "ymax") == 0){
            ymax = double_tmp;
          }
          else {printf("Error reading ymax.\n"); }

          fgets(file_line, 256, in_file);
          while(file_line[0] == '#'){
            fgets(file_line, 256, in_file);
          }
          sscanf(file_line, "%s %lf", key_name, &double_tmp);
          if(strcmp(key_name, "zmin") == 0){
            zmin = double_tmp;
          }
          else {printf("Error reading zmin.\n"); }

          fgets(file_line, 256, in_file);
          while(file_line[0] == '#'){
            fgets(file_line, 256, in_file);
          }
          sscanf(file_line, "%s %lf", key_name, &double_tmp);
          if(strcmp(key_name, "zmax") == 0){
            zmax = double_tmp;
          }
          else {printf("Error reading zmax.\n"); }
        }
        else if(strcmp(key_name, "Space_Grid") == 0){
          fgets(file_line, 256, in_file);
          while(file_line[0] == '#'){
            fgets(file_line, 256, in_file);
          }
          sscanf(file_line, "%s %d", key_name, &int_tmp);
          if(strcmp(key_name, "nx") == 0){
            nx = int_tmp;
          }
          else {printf("Error reading nx.\n"); }

          fgets(file_line, 256, in_file);
          while(file_line[0] == '#'){
            fgets(file_line, 256, in_file);
          }
          sscanf(file_line, "%s %d", key_name, &int_tmp);
          if(strcmp(key_name, "ny") == 0){
            ny = int_tmp;
          }
          else {printf("Error reading ny.\n"); }

          fgets(file_line, 256, in_file);
          while(file_line[0] == '#'){
            fgets(file_line, 256, in_file);
          }
          sscanf(file_line, "%s %d", key_name, &int_tmp);
          if(strcmp(key_name, "nz") == 0){
            nz = int_tmp;
          }
          else {printf("Error reading nz.\n"); }
        }
        else if(strcmp(key_name, "DiscreteOrdinates_Grid") == 0){
          fgets(file_line, 256, in_file);
          while(file_line[0] == '#'){
            fgets(file_line, 256, in_file);
          }
          sscanf(file_line, "%s %d", key_name, &int_tmp);
          if(strcmp(key_name, "num_dirsets_per_octant") == 0){
            num_dirsets_per_octant = int_tmp;
          }
          else {printf("Error reading num_dirsets_per_octant.\n"); }

          fgets(file_line, 256, in_file);
          while(file_line[0] == '#'){
            fgets(file_line, 256, in_file);
          }
          sscanf(file_line, "%s %d", key_name, &int_tmp);
          if(strcmp(key_name, "num_dirs_per_dirset") == 0){
            num_dirs_per_dirset = int_tmp;
          }
          else {printf("Error reading num_dir_per_dirset.\n"); }
        }
        else if(strcmp(key_name, "Energy_Grid") == 0){
          fgets(file_line, 256, in_file);
          while(file_line[0] == '#'){
            fgets(file_line, 256, in_file);
          }
          sscanf(file_line, "%s %d", key_name, &int_tmp);
          if(strcmp(key_name, "num_groupsets") == 0){
            num_groupsets = int_tmp;
          }
          else {printf("Error reading num_groupsets.\n"); }

          fgets(file_line, 256, in_file);
          while(file_line[0] == '#'){
            fgets(file_line, 256, in_file);
          }
          sscanf(file_line, "%s %d", key_name, &int_tmp);
          if(strcmp(key_name, "num_groups_per_groupset") == 0){
            num_groups_per_groupset = int_tmp;
          }
          else {printf("Error reading num_groupsets.\n"); }
        }
        else if(strcmp(key_name, "Boundary_Condition_Types") == 0){
          fgets(file_line, 256, in_file);
          while(file_line[0] == '#'){
            fgets(file_line, 256, in_file);
          }
          sscanf(file_line, "%s %d", key_name, &int_tmp);
          sprintf(test_key, "flux_left_type");
          if(strcmp(key_name, test_key) == 0){
            bndry_types[0] = int_tmp;
          }
          else {printf("Error reading %s.\n", test_key); }

          fgets(file_line, 256, in_file);
          while(file_line[0] == '#'){
            fgets(file_line, 256, in_file);
          }
          sscanf(file_line, "%s %d", key_name, &int_tmp);
          sprintf(test_key, "flux_right_type");
          if(strcmp(key_name, test_key) == 0){
            bndry_types[1] = int_tmp;
          }
          else {printf("Error reading %s.\n", test_key); }

          fgets(file_line, 256, in_file);
          while(file_line[0] == '#'){
            fgets(file_line, 256, in_file);
          }
          sscanf(file_line, "%s %d", key_name, &int_tmp);
          sprintf(test_key, "flux_front_type");
          if(strcmp(key_name, test_key) == 0){
            bndry_types[2] = int_tmp;
          }
          else {printf("Error reading %s.\n", test_key); }

          fgets(file_line, 256, in_file);
          while(file_line[0] == '#'){
            fgets(file_line, 256, in_file);
          }
          sscanf(file_line, "%s %d", key_name, &int_tmp);
          sprintf(test_key, "flux_back_type");
          if(strcmp(key_name, test_key) == 0){
            bndry_types[3] = int_tmp;
          }
          else {printf("Error reading %s.\n", test_key); }

          fgets(file_line, 256, in_file);
          while(file_line[0] == '#'){
            fgets(file_line, 256, in_file);
          }
          sscanf(file_line, "%s %d", key_name, &int_tmp);
          sprintf(test_key, "flux_lower_type");
          if(strcmp(key_name, test_key) == 0){
            bndry_types[4] = int_tmp;
          }
          else {printf("Error reading %s.\n", test_key); }

          fgets(file_line, 256, in_file);
          while(file_line[0] == '#'){
            fgets(file_line, 256, in_file);
          }
          sscanf(file_line, "%s %d", key_name, &int_tmp);
          sprintf(test_key, "flux_upper_type");
          if(strcmp(key_name, test_key) == 0){
            bndry_types[5] = int_tmp;
          }
          else {printf("Error reading %s.\n", test_key); }

        }
        else if(strcmp(key_name, "Boundary_Condition_Values") == 0){
          fgets(file_line, 256, in_file);
          while(file_line[0] == '#'){
            fgets(file_line, 256, in_file);
          }
          sscanf(file_line, "%s %lf", key_name, &double_tmp);
          sprintf(test_key, "flux_left_value");
          if(strcmp(key_name, test_key) == 0){
            bndry_values[0] = double_tmp;
          }

          fgets(file_line, 256, in_file);
          while(file_line[0] == '#'){
            fgets(file_line, 256, in_file);
          }
          sscanf(file_line, "%s %lf", key_name, &double_tmp);
          sprintf(test_key, "flux_right_value");
          if(strcmp(key_name, test_key) == 0){
            bndry_values[1] = double_tmp;
          }
          else {printf("Error reading %s.\n", test_key); }

          fgets(file_line, 256, in_file);
          while(file_line[0] == '#'){
            fgets(file_line, 256, in_file);
          }
          sscanf(file_line, "%s %lf", key_name, &double_tmp);
          sprintf(test_key, "flux_front_value");
          if(strcmp(key_name, test_key) == 0){
            bndry_values[2] = double_tmp;
          }
          else {printf("Error reading %s.\n", test_key); }

          fgets(file_line, 256, in_file);
          while(file_line[0] == '#'){
            fgets(file_line, 256, in_file);
          }
          sscanf(file_line, "%s %lf", key_name, &double_tmp);
          sprintf(test_key, "flux_back_value");
          if(strcmp(key_name, test_key) == 0){
            bndry_values[3] = double_tmp;
          }
          else {printf("Error reading %s.\n", test_key); }

          fgets(file_line, 256, in_file);
          while(file_line[0] == '#'){
            fgets(file_line, 256, in_file);
          }
          sscanf(file_line, "%s %lf", key_name, &double_tmp);
          sprintf(test_key, "flux_lower_value");
          if(strcmp(key_name, test_key) == 0){
            bndry_values[4] = double_tmp;
          }
          else {printf("Error reading %s.\n", test_key); }

          fgets(file_line, 256, in_file);
          while(file_line[0] == '#'){
            fgets(file_line, 256, in_file);
          }
          sscanf(file_line, "%s %lf", key_name, &double_tmp);
          sprintf(test_key, "flux_upper_value");
          if(strcmp(key_name, test_key) == 0){
            bndry_values[5] = double_tmp;
          }
          else {printf("Error reading %s.\n", test_key); }

        }
        else if(strcmp(key_name, "Source") == 0){
          fgets(file_line, 256, in_file);
          while(file_line[0] == '#'){
            fgets(file_line, 256, in_file);
          }
          sscanf(file_line, "%s %lf", key_name, &double_tmp);
          sprintf(test_key, "source_value");
          if(strcmp(key_name, test_key) == 0){
            source_value = double_tmp;
          }
          else {printf("Error reading source_value.\n"); }
        }
        else if(strcmp(key_name, "Sigma_Total") == 0){
          fgets(file_line, 256, in_file);
          while(file_line[0] == '#'){
            fgets(file_line, 256, in_file);
          }
          sscanf(file_line, "%s %lf", key_name, &double_tmp);
          sprintf(test_key, "sigma_total_value");
          if(strcmp(key_name, test_key) == 0){
            sigma_total_value = double_tmp;
          }
          else {printf("Error reading sigma_total_value.\n"); }
        }
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

  printf("Run_Name\n");
  printf("   run_name = %s\n", run_name);
  printf("\n");
  printf("Processors\n");
  printf("   npx = %i   npy = %i   npz = %i\n",
          npx, npy, npz);
  printf("   ncalls = %i\n", ncalls);
  printf("\n");
  printf("Domain\n");
  printf("   xmin = %e  xmax = %e\n", xmin,
          xmax);
  printf("   ymin = %e  ymax = %e\n", ymin,
          ymax);
  printf("   zmin = %e  zmax = %e\n", zmin,
          zmax);
  printf("\n");
  printf("Space_Grid\n");
  printf("   nx = %i  ny = %i  nz = %i\n",
          nx, ny, nz);
  printf("\n");
  printf("DiscreteOrdinates_Grid\n");
  printf("   num_dirsets_per_octant =  %i\n", num_dirsets_per_octant);
  printf("   num_dirs_per_dirset =      %i\n", num_dirs_per_dirset);
  printf("\n");
  printf("Energy_Grid\n");
  printf("   num_groupsets =           %i\n", num_groupsets);
  printf("   num_groups_per_groupset = %i\n", num_groups_per_groupset);
  printf("Boundary_Conditions\n");
  for(int j = 0; j < 6; j++){
    if( (bndry_types[j]) == 0){
      printf("      Dirichlet  Value = %e\n",
              bndry_values[j]);
    }
    else if( (bndry_types[j]) == 1){
      printf("      Reflecting Value = %e\n",
              bndry_values[j]);
    }
  }
  printf("\n");
  printf("Source\n");
  printf("   source_value = %e\n", source_value);
  printf("\n");
  printf("Sigma_Total\n");
  printf("   sigma_total_value = %e\n",
          sigma_total_value);

}

