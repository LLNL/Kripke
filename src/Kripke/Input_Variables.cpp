/*
 * NOTICE
 *
 * This work was produced at the Lawrence Livermore National Laboratory (LLNL)
 * under contract no. DE-AC-52-07NA27344 (Contract 44) between the U.S.
 * Department of Energy (DOE) and Lawrence Livermore National Security, LLC
 * (LLNS) for the operation of LLNL. The rights of the Federal Government are
 * reserved under Contract 44.
 *
 * DISCLAIMER
 *
 * This work was prepared as an account of work sponsored by an agency of the
 * United States Government. Neither the United States Government nor Lawrence
 * Livermore National Security, LLC nor any of their employees, makes any
 * warranty, express or implied, or assumes any liability or responsibility
 * for the accuracy, completeness, or usefulness of any information, apparatus,
 * product, or process disclosed, or represents that its use would not infringe
 * privately-owned rights. Reference herein to any specific commercial products,
 * process, or service by trade name, trademark, manufacturer or otherwise does
 * not necessarily constitute or imply its endorsement, recommendation, or
 * favoring by the United States Government or Lawrence Livermore National
 * Security, LLC. The views and opinions of authors expressed herein do not
 * necessarily state or reflect those of the United States Government or
 * Lawrence Livermore National Security, LLC, and shall not be used for
 * advertising or product endorsement purposes.
 *
 * NOTIFICATION OF COMMERCIAL USE
 *
 * Commercialization of this product is prohibited without notifying the
 * Department of Energy (DOE) or Lawrence Livermore National Security.
 */

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
  num_groupsets(2),
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
