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

#include<Kripke/Kernel.h>
#include<Kripke/Grid.h>
#include<Kripke/SubTVec.h>

#include<Domain/Layout.h>
#include<Domain/Forall.h>


#include<Kripke/Kernel/Kernel_3d_GDZ.h>
#include<Kripke/Kernel/Kernel_3d_DGZ.h>
#include<Kripke/Kernel/Kernel_3d_ZDG.h>
#include<Kripke/Kernel/Kernel_3d_DZG.h>
#include<Kripke/Kernel/Kernel_3d_ZGD.h>
#include<Kripke/Kernel/Kernel_3d_GZD.h>

#include<Kripke/Kernel/VariablePolicy.h>

/**
 * Factory to create a kernel object for the specified nesting
 */
Kernel *createKernel(Nesting_Order nest, int num_dims){
  if(num_dims == 3){
    switch(nest){
    case NEST_GDZ:
      return new Kernel_3d_GDZ();
    case NEST_DGZ:
      return new Kernel_3d_DGZ();
    case NEST_ZDG:
      return new Kernel_3d_ZDG();
    case NEST_DZG:
      return new Kernel_3d_DZG();
    case NEST_ZGD:
      return new Kernel_3d_ZGD();
    case NEST_GZD:
      return new Kernel_3d_GZD();
    }
  }

  MPI_Abort(MPI_COMM_WORLD, 1);
  return NULL;
}


Kernel::Kernel(Nesting_Order nest) :
  nesting_order(nest)
{}

Kernel::~Kernel(){
}




#include<Kripke/Kernel/LTimesPolicy.h>
void Kernel::LTimes(Grid_Data *grid_data) {

  getNestType(nesting_order, nest_type, [=](){
    typedef VariablePolicy<nest_type> VPOL;
    typedef VariableView<VPOL> VIEW;
 
    // Outer parameters
    int num_moments = grid_data->total_num_moments;

    // Zero Phi
    for(int ds = 0;ds < grid_data->num_zone_sets;++ ds){
      grid_data->phi[ds]->clear(0.0);
    }

    // Loop over Subdomains
    int num_subdomains = grid_data->subdomains.size();
    for (int sdom_id = 0; sdom_id < num_subdomains; ++ sdom_id){
      Subdomain &sdom = grid_data->subdomains[sdom_id];
    
      // Get dimensioning
      int num_zones = sdom.num_zones;
      int num_groups = sdom.phi->groups;
      int num_local_groups = sdom.num_groups;
      int group0 = sdom.group0;
      int num_local_directions = sdom.num_directions;
                   
      // Get pointers    
      typename VIEW::View3d_Psi psi(sdom.psi->ptr(), num_local_directions, num_local_groups, num_zones);
      typename VIEW::View3d_Phi phi(sdom.phi->ptr(), num_moments, num_groups, num_zones);
      typename VIEW::View2d_Ell ell(sdom.ell->ptr(), num_local_directions, num_moments);
            
      forall4<LTimesPolicy<nest_type> >(
        num_moments, num_local_directions, num_local_groups, num_zones, 
        [&](int nm, int d, int g, int z){
 
          phi(nm, g+group0, z) += ell(d,nm) * psi(d,g,z);
 
        });
    }
  });
}

#include<Kripke/Kernel/LPlusTimesPolicy.h>
void Kernel::LPlusTimes(Grid_Data *grid_data) {

  getNestType(nesting_order, nest_type, [&](){
    typedef VariablePolicy<nest_type> VPOL;
    typedef VariableView<VPOL> VIEW;
 
    // Outer parameters
    int num_moments = grid_data->total_num_moments;

    // Loop over Subdomains
    int num_subdomains = grid_data->subdomains.size();
    for (int sdom_id = 0; sdom_id < num_subdomains; ++ sdom_id){
      Subdomain &sdom = grid_data->subdomains[sdom_id];
      // Zero RHS
      sdom.rhs->clear(0.0);
    }
    for (int sdom_id = 0; sdom_id < num_subdomains; ++ sdom_id){
      Subdomain &sdom = grid_data->subdomains[sdom_id];

      // Get dimensioning
      int num_zones = sdom.num_zones;
      int num_local_groups = sdom.num_groups;
      int num_groups = sdom.phi_out->groups;
      int group0 = sdom.group0;
      int num_local_directions = sdom.num_directions;
      
      // Get pointers
      VIEW::View3d_Psi     rhs(sdom.rhs->ptr(), num_local_directions, num_local_groups, num_zones);
      VIEW::View3d_Phi     phi_out(sdom.phi_out->ptr(), num_moments, num_groups, num_zones);
      VIEW::View2d_EllPlus ell_plus(sdom.ell_plus->ptr(), num_local_directions, num_moments);

      forall4<LPlusTimesPolicy<nest_type> >(
        num_moments, num_local_directions, num_local_groups, num_zones, 
        [&](int nm, int d, int g, int z){
 
          rhs(d,g,z) += ell_plus(d,nm) * phi_out(nm,g+group0,z);
 
        });
    }
  });
}



