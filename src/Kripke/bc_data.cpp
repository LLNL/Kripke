/*--------------------------------------------------------------------------
 * Utility functions for BC_Data structure
 *--------------------------------------------------------------------------*/

#include "transport_headers.h"

/*--------------------------------------------------------------------------
 * NewBCData : Creates a new BC_Data and allocates memory for
 *             its data.
 *--------------------------------------------------------------------------*/

BC_Data *NewBCData(Grid_Data *grid_data)
/*--------------------------------------------------------------------------
 * grid_data : Grid information structure
 *--------------------------------------------------------------------------*/
{
  BC_Data *bc_data;

  int     *nzones         = grid_data->nzones;

  int ndir;
  int dirmin, dirmax;
  int n1, n2;
  int face;
  int g;

  NEW(bc_data, 1, BC_Data *);

  return(bc_data);
}

/*--------------------------------------------------------------------------
 * InitBCData : Initializes an BC_Data structure given a set
 *              of boundary conditions.  It is assumed that boundary
 *              conditions are constant on a face and either of
 *              Dirichlet or Reflecting type.
 *--------------------------------------------------------------------------*/

void InitBCData(int *types, double *vals, Grid_Data *grid_data,
                BC_Data *bc_data)
/*--------------------------------------------------------------------------
 * types     : Array of integers indicating the boundary condition type
 *             on each face for each energy group.
 *                            Type = 0  => Dirichlet
 *                            Type = 1  => Specular
 * vals      : Array of doubles with boundary data
 * grid_data : Grid information
 * bc_data   : The BC_Data structure of boundary condition
 *             information to be returned
 *--------------------------------------------------------------------------*/
{
  int  *nzones       = grid_data->nzones;
  int i1, i2;
  int ndir;
  int dirmin, dirmax;
  int n1, n2;
  int face;
  int index;
  int bc_type;

  for(ndir = 0; ndir < 3; ndir++){
    if( ((ndir+1)%3) > ((ndir+2)%3) ){
      dirmin = (ndir+2)%3;
      dirmax = (ndir+1)%3;
    }
    else {
      dirmin = (ndir+1)%3;
      dirmax = (ndir+2)%3;
    }
    n1 = nzones[dirmin];
    n2 = nzones[dirmax];

    for(face = 0; face < 2; face++){
      if( (grid_data->mynbr[ndir][face]) == -1){
        bc_data->bc_types[ndir*2 + face] = types[ndir*2 + face];
        bc_type = types[ndir*2 + face];

        if( (bc_type == 0) || (bc_type == 1) ){
          /* Dirichlet or Neumann condition */
          bc_data->bc_values[ndir*2 + face] = vals[ndir*2 + face];
        }
        else {
          /* Illegal type */
          error_exit(1);
        }
      }
    }
  }

}

/*--------------------------------------------------------------------------
 * FreeBCData : Frees all memory associated with an BC_Data
 *                    structure
 *--------------------------------------------------------------------------*/

void FreeBCData(BC_Data *bc_data)
/*--------------------------------------------------------------------------
 * bc_data : The BC_Data structure to be deallocated
 *--------------------------------------------------------------------------*/
{
  FREE(bc_data);
}
