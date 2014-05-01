/*--------------------------------------------------------------------------
 * Initialization and evaluation functions for Sigma_Tot.
 *--------------------------------------------------------------------------*/

#include "transport_headers.h"

/*--------------------------------------------------------------------------
 * NewSigmaTot : Creates a new Sigma_Tot structure and
 *               allocates memory for it.
 *--------------------------------------------------------------------------*/

Sigma_Tot *NewSigmaTot(double param0)
/*--------------------------------------------------------------------------
 * type       : The type of cross section requested
 *              type = 0  => constant
 * param0     : For type = 0, the constant value.
 *--------------------------------------------------------------------------*/
{
  Sigma_Tot          *sigma_tot;

  NEW(sigma_tot, 1, Sigma_Tot *);
  sigma_tot->data = param0;

  return(sigma_tot);
}


/*--------------------------------------------------------------------------
 * FreeSigmaTot : Frees a Sigma_Tot structure
 *--------------------------------------------------------------------------*/

void FreeSigmaTot(Sigma_Tot *sigma_tot)
/*--------------------------------------------------------------------------
 * sigma_tot : The Sigma_Tot structure to be deallocated.
 *--------------------------------------------------------------------------*/
{
  FREE(sigma_tot);
}

/*--------------------------------------------------------------------------
 * EvalSigmaTot : Evaluate the sigma_tot function for each spatial zone
 *--------------------------------------------------------------------------*/

void EvalSigmaTot(Sigma_Tot *sigma_tot, Data_Vector *vector)
/*--------------------------------------------------------------------------
 * sigma_tot    : The Sigma_Tot structure that specifies the
 *                cross section function
 * group        : The energy group where this function should be evaluated
 * vector       : The vector of calculated values to be returned
 *--------------------------------------------------------------------------*/
{
  int   *nzones = vector->grid_data->nzones;
  int num_zones = nzones[0]*nzones[1]*nzones[2];
  int zone;

  /* Constant cross section case */
  double value = sigma_tot->data;

  for(zone=0; zone<num_zones; zone++){
    vector->data[zone] = value;
  }

}
