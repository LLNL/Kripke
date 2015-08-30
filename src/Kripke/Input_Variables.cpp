#include<Kripke/Input_Variables.h>

/**
* Setup the default input choices
*/
Input_Variables::Input_Variables() : 
  run_name("kripke"),
  nx(16), ny(16), nz(16),
  num_directions(96),
  num_groups(32),
  legendre_order(4),
  quad_num_polar(0),
  quad_num_azimuthal(0),
 
  nesting(NEST_DGZ),
 
  npx(1), npy(1), npz(1),
  num_dirsets(8),
  num_groupsets(1),
  layout_pattern(0),
  
  niter(10),
  parallel_method(PMETHOD_SWEEP)
{
  num_zonesets_dim[0] = 1; 
  num_zonesets_dim[1] = 1;
  num_zonesets_dim[2] = 1;

  sigt[0] = 0.1;  
  sigt[1] = 0.0001;
  sigt[2] = 0.1;
  
  sigs[0] = 0.05;  
  sigs[1] = 0.00005;
  sigs[2] = 0.05; 
}

/**
 *  Checks validity of inputs, returns 'true' on error.
 */
bool Input_Variables::checkValues(void) const{
  if(num_zonesets_dim[0] <= 0 || num_zonesets_dim[1] <= 0 || num_zonesets_dim[2] <= 0){
    return true;
  }
  
  if(layout_pattern < 0 || layout_pattern > 1){
    return true;
  }
  
  if(nesting < 0){
    return true;
  }
  
  return false;
}
