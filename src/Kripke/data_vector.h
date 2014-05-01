/*--------------------------------------------------------------------------
 * Header file for the Data_Vector data structures
 *--------------------------------------------------------------------------*/

#ifndef included_data_vector
#define included_data_vector

#include "grid.h"

/*--------------------------------------------------------------------------
 * Define the Data_Vector structure.
 *
 * data      : Vector of length equal to the number of cells.
 * grid_data : Grid information structure
 *--------------------------------------------------------------------------*/

typedef struct {
  Grid_Data *grid_data;

  double    *data;

}   Data_Vector;

#endif
