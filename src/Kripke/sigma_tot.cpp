
#include<Kripke/user_data.h>
#include<vector>

/*--------------------------------------------------------------------------
 * EvalSigmaTot : Evaluate the sigma_tot function for each spatial zone
 *--------------------------------------------------------------------------*/

void EvalSigmaTot(User_Data *user_data, std::vector<double> &vector)
/*--------------------------------------------------------------------------
 * sigma_tot    : The Sigma_Tot structure that specifies the
 *                cross section function
 * group        : The energy group where this function should be evaluated
 * vector       : The vector of calculated values to be returned
 *--------------------------------------------------------------------------*/
{
  int   *nzones = user_data->grid_data->nzones;
  int num_zones = nzones[0]*nzones[1]*nzones[2];
  int zone;

  /* Constant cross section case */
  double value = user_data->sigma_tot;

  vector.resize(num_zones);
  for(zone=0; zone<num_zones; zone++){
    vector[zone] = value;
  }

}
