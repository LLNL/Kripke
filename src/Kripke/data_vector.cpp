/*--------------------------------------------------------------------------
 * Utility functions for the Data_Vector structure.
 *--------------------------------------------------------------------------*/

#include "transport_headers.h"


/*--------------------------------------------------------------------------
 * NewDataVector : Creates a new Data_Vector and allocates
 *                       memory for its data.
 *--------------------------------------------------------------------------*/

Data_Vector *NewDataVector(Grid_Data *grid_data)
/*--------------------------------------------------------------------------
 * grid_data : Grid information structure
 *--------------------------------------------------------------------------*/
{
  int length = (grid_data->nzones[0]) * (grid_data->nzones[1]) *
               (grid_data->nzones[2]);
  Data_Vector *data_vector;

  NEW(data_vector, 1, Data_Vector *);
  NEW(data_vector->data, length, double *);

  data_vector->grid_data = grid_data;

  return(data_vector);
}

/*--------------------------------------------------------------------------
 * FreeDataVector : Frees all memory associated with a Data_Vector
 *--------------------------------------------------------------------------*/

void FreeDataVector(Data_Vector *data_vector)
/*--------------------------------------------------------------------------
 * data_vector    : The Data_Vector to be deallocated.
 *--------------------------------------------------------------------------*/
{
  FREE(data_vector->data);
  FREE(data_vector);
}
