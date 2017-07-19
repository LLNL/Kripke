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

#include <Kripke/Grid.h>
#include <Kripke/SubTVec.h>
#include <Kripke/InputVariables.h>

Subdomain::Subdomain()
{

}
Subdomain::~Subdomain(){

}


/**
  Setup subdomain and allocate data
*/
void Subdomain::setup(int sdom_id, InputVariables *input_vars, int , int ds, int ,
    std::vector<Kripke::QuadraturePoint> &direction_list, Layout *layout)
{

  int num_directions = input_vars->num_directions / input_vars->num_dirsets;
  int direction0 = ds * num_directions;


  // Setup neighbor data
  int dirs[3] = {
      direction_list[direction0].id,
      direction_list[direction0].jd,
      direction_list[direction0].kd};
  for(int dim = 0;dim < 3;++ dim){
    downwind[dim] = layout->getNeighbor(sdom_id, dim, dirs[dim]);
    upwind[dim] = layout->getNeighbor(sdom_id, dim, -1 * dirs[dim]);
  }


}






