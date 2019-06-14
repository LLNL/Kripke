//
// Copyright (c) 2014-19, Lawrence Livermore National Security, LLC
// and Kripke project contributors. See the COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
//

#include <Kripke/InputVariables.h>

#include <Kripke/Core/Comm.h>

using namespace Kripke;

/**
* Setup the default input choices
*/
InputVariables::InputVariables() : 
  nx(16), ny(16), nz(16),
  num_directions(96),
  num_groups(32),
  legendre_order(4),
  quad_num_polar(0),
  quad_num_azimuthal(0),
 
  al_v(ArchLayoutV{KRIPKE_ARCHV_DEFAULT, KRIPKE_LAYOUTV_DEFAULT}),
 
  npx(1), npy(1), npz(1),
  num_dirsets(8),
  num_groupsets(2),
  
  niter(10),
  parallel_method(PMETHOD_SWEEP),
  num_material_subsamples(4),
  run_name("kripke")
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
bool InputVariables::checkValues(void) const{
  // make sure any output only goes to root

  Kripke::Core::Comm comm;
  int rank = comm.rank();

  if(num_zonesets_dim[0] <= 0 || num_zonesets_dim[1] <= 0 || num_zonesets_dim[2] <= 0){
    if(!rank)
      printf("Number of zone-sets in each dim need to be greater than or equal to 1\n");
    return true;
  }
  
  if(num_groups < 1){
    if(!rank)
      printf("Number of groups (%d) needs to be at least 1\n", num_groups);
    return true;
  }
  
  if(num_groups % num_groupsets){
    if(!rank)
      printf("Number of groups (%d) must be evenly divided by number of groupsets (%d)\n",
        num_groups, num_groupsets);
    return true;
  }
  
  if(num_directions < 8){
    if(!rank)
      printf("Number of directions (%d) needs to be at least 8\n", num_directions);
    return true;
  }
  
  if(num_dirsets % 8 && num_dirsets < 8){
    if(!rank)
      printf("Number of direction sets (%d) must be a multiple of 8\n", num_dirsets);
    return true;
  }
  
  if(num_directions % num_dirsets){
    if(!rank)
      printf("Number of directions (%d) must be evenly divided by number of directionsets(%d)\n",
        num_directions, num_dirsets);
    return true;
  }
  
  if(legendre_order < 0){
    if(!rank)
      printf("Legendre scattering order (%d) must be >= 0\n", legendre_order);
    return true;
  }
  
  if(niter < 1){
    if(!rank)
      printf("You must run at least one iteration (%d)\n", niter);
    return true;
  }
  
  return false;
}
