/*--------------------------------------------------------------------------
 * Input-reading utilities
 *--------------------------------------------------------------------------*/

#include "input_variables.h"
#include "transport_headers.h"

#include <stdio.h>

/*--------------------------------------------------------------------------
 * Utility routines for the Input_Variables structure.
 *--------------------------------------------------------------------------*/

Input_Variables *ReadInput(FILE *in_file)
/*--------------------------------------------------------------------------
 * in_file: Input file
 *--------------------------------------------------------------------------*/
{
  Input_Variables *input_variables;

  char file_line[256];
  char key_name[256];
  char test_key[256];
  char str_tmp[256];

  int n, i, j;

  int int_tmp;
  double double_tmp;

  NEW(input_variables, 1, Input_Variables *);

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
          strcpy(input_variables->run_name, str_tmp);
        }
      }
      else if(strcmp(key_name, "Processors") == 0){
        fgets(file_line, 256, in_file);
        while(file_line[0] == '#'){
          fgets(file_line, 256, in_file);
        }
        sscanf(file_line, "%s %d", key_name, &int_tmp);
        if(strcmp(key_name, "npx") == 0){
          input_variables->npx = int_tmp;
        }
        else {printf("Error reading npx.\n"); }

        fgets(file_line, 256, in_file);
        while(file_line[0] == '#'){
          fgets(file_line, 256, in_file);
        }
        sscanf(file_line, "%s %d", key_name, &int_tmp);
        if(strcmp(key_name, "npy") == 0){
          input_variables->npy = int_tmp;
        }
        else {printf("Error reading npy.\n"); }

        fgets(file_line, 256, in_file);
        while(file_line[0] == '#'){
          fgets(file_line, 256, in_file);
        }
        sscanf(file_line, "%s %d", key_name, &int_tmp);
        if(strcmp(key_name, "npz") == 0){
          input_variables->npz = int_tmp;
        }
        else {printf("Error reading npz.\n"); }

        fgets(file_line, 256, in_file);
        while(file_line[0] == '#'){
          fgets(file_line, 256, in_file);
        }
        sscanf(file_line, "%s %d", key_name, &int_tmp);
        if(strcmp(key_name, "nlevels_kba") == 0){
          input_variables->nlevels_kba = int_tmp;
        }
        else {printf("Error reading nlevels_kba.\n"); }

        fgets(file_line, 256, in_file);
        while(file_line[0] == '#'){
          fgets(file_line, 256, in_file);
        }
        sscanf(file_line, "%s %d", key_name, &int_tmp);
        if(strcmp(key_name, "ncalls") == 0){
          input_variables->ncalls = int_tmp;
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
          input_variables->xmin = double_tmp;
        }
        else {printf("Error reading xmin.\n"); }

        fgets(file_line, 256, in_file);
        while(file_line[0] == '#'){
          fgets(file_line, 256, in_file);
        }
        sscanf(file_line, "%s %lf", key_name, &double_tmp);
        if(strcmp(key_name, "xmax") == 0){
          input_variables->xmax = double_tmp;
        }
        else {printf("Error reading xmax.\n"); }

        fgets(file_line, 256, in_file);
        while(file_line[0] == '#'){
          fgets(file_line, 256, in_file);
        }
        sscanf(file_line, "%s %lf", key_name, &double_tmp);
        if(strcmp(key_name, "ymin") == 0){
          input_variables->ymin = double_tmp;
        }
        else {printf("Error reading ymin.\n"); }

        fgets(file_line, 256, in_file);
        while(file_line[0] == '#'){
          fgets(file_line, 256, in_file);
        }
        sscanf(file_line, "%s %lf", key_name, &double_tmp);
        if(strcmp(key_name, "ymax") == 0){
          input_variables->ymax = double_tmp;
        }
        else {printf("Error reading ymax.\n"); }

        fgets(file_line, 256, in_file);
        while(file_line[0] == '#'){
          fgets(file_line, 256, in_file);
        }
        sscanf(file_line, "%s %lf", key_name, &double_tmp);
        if(strcmp(key_name, "zmin") == 0){
          input_variables->zmin = double_tmp;
        }
        else {printf("Error reading zmin.\n"); }

        fgets(file_line, 256, in_file);
        while(file_line[0] == '#'){
          fgets(file_line, 256, in_file);
        }
        sscanf(file_line, "%s %lf", key_name, &double_tmp);
        if(strcmp(key_name, "zmax") == 0){
          input_variables->zmax = double_tmp;
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
          input_variables->nx = int_tmp;
        }
        else {printf("Error reading nx.\n"); }

        fgets(file_line, 256, in_file);
        while(file_line[0] == '#'){
          fgets(file_line, 256, in_file);
        }
        sscanf(file_line, "%s %d", key_name, &int_tmp);
        if(strcmp(key_name, "ny") == 0){
          input_variables->ny = int_tmp;
        }
        else {printf("Error reading ny.\n"); }

        fgets(file_line, 256, in_file);
        while(file_line[0] == '#'){
          fgets(file_line, 256, in_file);
        }
        sscanf(file_line, "%s %d", key_name, &int_tmp);
        if(strcmp(key_name, "nz") == 0){
          input_variables->nz = int_tmp;
        }
        else {printf("Error reading nz.\n"); }
      }
      else if(strcmp(key_name, "DiscreteOrdinates_Grid") == 0){
        fgets(file_line, 256, in_file);
        while(file_line[0] == '#'){
          fgets(file_line, 256, in_file);
        }
        sscanf(file_line, "%s %d", key_name, &int_tmp);
        if(strcmp(key_name, "num_directions_per_octant") == 0){
          input_variables->num_directions_per_octant = int_tmp;
        }
        else {printf("Error reading num_directions_per_octant.\n"); }
      }
      else if(strcmp(key_name, "Boundary_Condition_Types") == 0){
        NEW(input_variables->bndry_types, 6, int *);
        fgets(file_line, 256, in_file);
        while(file_line[0] == '#'){
          fgets(file_line, 256, in_file);
        }
        sscanf(file_line, "%s %d", key_name, &int_tmp);
        sprintf(test_key, "flux_left_type");
        if(strcmp(key_name, test_key) == 0){
          input_variables->bndry_types[0] = int_tmp;
        }
        else {printf("Error reading %s.\n", test_key); }

        fgets(file_line, 256, in_file);
        while(file_line[0] == '#'){
          fgets(file_line, 256, in_file);
        }
        sscanf(file_line, "%s %d", key_name, &int_tmp);
        sprintf(test_key, "flux_right_type");
        if(strcmp(key_name, test_key) == 0){
          input_variables->bndry_types[1] = int_tmp;
        }
        else {printf("Error reading %s.\n", test_key); }

        fgets(file_line, 256, in_file);
        while(file_line[0] == '#'){
          fgets(file_line, 256, in_file);
        }
        sscanf(file_line, "%s %d", key_name, &int_tmp);
        sprintf(test_key, "flux_front_type");
        if(strcmp(key_name, test_key) == 0){
          input_variables->bndry_types[2] = int_tmp;
        }
        else {printf("Error reading %s.\n", test_key); }

        fgets(file_line, 256, in_file);
        while(file_line[0] == '#'){
          fgets(file_line, 256, in_file);
        }
        sscanf(file_line, "%s %d", key_name, &int_tmp);
        sprintf(test_key, "flux_back_type");
        if(strcmp(key_name, test_key) == 0){
          input_variables->bndry_types[3] = int_tmp;
        }
        else {printf("Error reading %s.\n", test_key); }

        fgets(file_line, 256, in_file);
        while(file_line[0] == '#'){
          fgets(file_line, 256, in_file);
        }
        sscanf(file_line, "%s %d", key_name, &int_tmp);
        sprintf(test_key, "flux_lower_type");
        if(strcmp(key_name, test_key) == 0){
          input_variables->bndry_types[4] = int_tmp;
        }
        else {printf("Error reading %s.\n", test_key); }

        fgets(file_line, 256, in_file);
        while(file_line[0] == '#'){
          fgets(file_line, 256, in_file);
        }
        sscanf(file_line, "%s %d", key_name, &int_tmp);
        sprintf(test_key, "flux_upper_type");
        if(strcmp(key_name, test_key) == 0){
          input_variables->bndry_types[5] = int_tmp;
        }
        else {printf("Error reading %s.\n", test_key); }

      }
      else if(strcmp(key_name, "Boundary_Condition_Values") == 0){
        NEW(input_variables->bndry_values, 6, double *);
        fgets(file_line, 256, in_file);
        while(file_line[0] == '#'){
          fgets(file_line, 256, in_file);
        }
        sscanf(file_line, "%s %lf", key_name, &double_tmp);
        sprintf(test_key, "flux_left_value");
        if(strcmp(key_name, test_key) == 0){
          input_variables->bndry_values[0] = double_tmp;
        }

        fgets(file_line, 256, in_file);
        while(file_line[0] == '#'){
          fgets(file_line, 256, in_file);
        }
        sscanf(file_line, "%s %lf", key_name, &double_tmp);
        sprintf(test_key, "flux_right_value");
        if(strcmp(key_name, test_key) == 0){
          input_variables->bndry_values[1] = double_tmp;
        }
        else {printf("Error reading %s.\n", test_key); }

        fgets(file_line, 256, in_file);
        while(file_line[0] == '#'){
          fgets(file_line, 256, in_file);
        }
        sscanf(file_line, "%s %lf", key_name, &double_tmp);
        sprintf(test_key, "flux_front_value");
        if(strcmp(key_name, test_key) == 0){
          input_variables->bndry_values[2] = double_tmp;
        }
        else {printf("Error reading %s.\n", test_key); }

        fgets(file_line, 256, in_file);
        while(file_line[0] == '#'){
          fgets(file_line, 256, in_file);
        }
        sscanf(file_line, "%s %lf", key_name, &double_tmp);
        sprintf(test_key, "flux_back_value");
        if(strcmp(key_name, test_key) == 0){
          input_variables->bndry_values[3] = double_tmp;
        }
        else {printf("Error reading %s.\n", test_key); }

        fgets(file_line, 256, in_file);
        while(file_line[0] == '#'){
          fgets(file_line, 256, in_file);
        }
        sscanf(file_line, "%s %lf", key_name, &double_tmp);
        sprintf(test_key, "flux_lower_value");
        if(strcmp(key_name, test_key) == 0){
          input_variables->bndry_values[4] = double_tmp;
        }
        else {printf("Error reading %s.\n", test_key); }

        fgets(file_line, 256, in_file);
        while(file_line[0] == '#'){
          fgets(file_line, 256, in_file);
        }
        sscanf(file_line, "%s %lf", key_name, &double_tmp);
        sprintf(test_key, "flux_upper_value");
        if(strcmp(key_name, test_key) == 0){
          input_variables->bndry_values[5] = double_tmp;
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
          input_variables->source_value = double_tmp;
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
          input_variables->sigma_total_value = double_tmp;
        }
        else {printf("Error reading sigma_total_value.\n"); }
      }
    }
  }

  return(input_variables);
}

/*--------------------------------------------------------------------------
 * PrintInputVariables : Prints values in an Input_Variables structure
 *--------------------------------------------------------------------------*/

void PrintInputVariables(Input_Variables *input_variables, FILE *out_file)
/*--------------------------------------------------------------------------
 * input_variables : The Input_Variables structure to be printed.
 * out_file        : The file to which input information will be printed.
 *--------------------------------------------------------------------------*/
{
  int j;

  fprintf(out_file, "Run_Name\n");
  fprintf(out_file, "   run_name = %s\n", input_variables->run_name);
  fprintf(out_file, "\n");
  fprintf(out_file, "Processors\n");
  fprintf(out_file, "   npx = %i   npy = %i   npz = %i   nlevels_kba = %i\n",
          input_variables->npx, input_variables->npy, input_variables->npz,
	  input_variables->nlevels_kba);
  fprintf(out_file, "   ncalls = %i\n", input_variables->ncalls);
  fprintf(out_file, "\n");
  fprintf(out_file, "Domain\n");
  fprintf(out_file, "   xmin = %e  xmax = %e\n", input_variables->xmin,
          input_variables->xmax);
  fprintf(out_file, "   ymin = %e  ymax = %e\n", input_variables->ymin,
          input_variables->ymax);
  fprintf(out_file, "   zmin = %e  zmax = %e\n", input_variables->zmin,
          input_variables->zmax);
  fprintf(out_file, "\n");
  fprintf(out_file, "Space_Grid\n");
  fprintf(out_file, "   nx = %i  ny = %i  nz = %i\n",
          input_variables->nx, input_variables->ny, input_variables->nz);
  fprintf(out_file, "\n");
  fprintf(out_file, "DiscreteOrdinates_Grid\n");
  fprintf(out_file, "   num_directions_per_octant = %i\n",
          input_variables->num_directions_per_octant);
  fprintf(out_file, "\n");
  fprintf(out_file, "Boundary_Conditions\n");
  for(j = 0; j < 6; j++){
    if( (input_variables->bndry_types[j]) == 0){
      fprintf(out_file, "      Dirichlet  Value = %e\n",
              input_variables->bndry_values[j]);
    }
    else if( (input_variables->bndry_types[j]) == 1){
      fprintf(out_file, "      Reflecting Value = %e\n",
              input_variables->bndry_values[j]);
    }
  }
  fprintf(out_file, "\n");
  fprintf(out_file, "Source\n");
  fprintf(out_file, "   source_value = %e\n", input_variables->source_value);
  fprintf(out_file, "\n");
  fprintf(out_file, "Sigma_Total\n");
  fprintf(out_file, "   sigma_total_value = %e\n",
          input_variables->sigma_total_value);

}


/*--------------------------------------------------------------------------
 * FreeInputVariables : Frees an Input_Variables structure
 *--------------------------------------------------------------------------*/

void FreeInputVariables(Input_Variables *input_variables)
/*--------------------------------------------------------------------------
 * input_variables : The Input_Variables structure to be deallocated.
 *--------------------------------------------------------------------------*/
{
  FREE(input_variables->bndry_types);
  FREE(input_variables->bndry_values);
  FREE(input_variables);
}
