#include<Kripke/Kernel.h>
#include<Kripke/Comm.h>
#include<Kripke/Grid.h>
#include<Kripke/User_Data.h>
#include<Kripke/SubTVec.h>

#include<Kripke/Kernel/Kernel_3d_GDZ.h>
#include<Kripke/Kernel/Kernel_3d_DGZ.h>
#include<Kripke/Kernel/Kernel_3d_ZDG.h>
#include<Kripke/Kernel/Kernel_3d_DZG.h>
#include<Kripke/Kernel/Kernel_3d_ZGD.h>
#include<Kripke/Kernel/Kernel_3d_GZD.h>

Kernel::Kernel(){

}

Kernel::~Kernel(){

}

/**
 * Allocates kernel-specific storage for the User_Data object.
 */
void Kernel::allocateStorage(User_Data *user_data){
  Grid_Data *grid_data = user_data->grid_data;
  Nesting_Order nest = nestingPsi();
  for(int gs = 0;gs < grid_data->gd_sets.size();++ gs){
    for(int ds = 0;ds < grid_data->gd_sets[gs].size();++ ds){
      Group_Dir_Set &gdset = grid_data->gd_sets[gs][ds];
      gdset.allocate(grid_data, nest);
    }
  }

  // Allocate moments variables
  int num_moments = grid_data->num_moments;
  int num_dims = 3;
  int total_dirs = user_data->directions.size();
  int num_zones = grid_data->num_zones;
  int num_groups = user_data->num_group_sets * user_data->num_groups_per_set;

  // Allocate L, L+ table nm indicies which make threading easier
  grid_data->nm_table.clear();
  for (int n = 0; n < num_moments; n++) {
    for (int m = -n; m <= n; m++) {
      grid_data->nm_table.push_back(n);
    }
  }

  int total_moments = grid_data->nm_table.size();
  grid_data->phi = new SubTVec(nestingPhi(), num_groups, total_moments, num_zones);
  grid_data->phi_out = new SubTVec(nestingPhi(), num_groups, total_moments, num_zones);

  if(nest == NEST_GDZ || nest == NEST_DZG || nest == NEST_DGZ){
    grid_data->ell = new SubTVec(NEST_ZGD, total_moments, total_dirs, 1);
    grid_data->ell_plus = new SubTVec(NEST_ZDG, total_moments, total_dirs, 1);
  }
  else{
    grid_data->ell = new SubTVec(NEST_ZDG, total_moments, total_dirs, 1);
    grid_data->ell_plus = new SubTVec(NEST_ZDG, total_moments, total_dirs, 1);
  }

  // allocate sigt  1xGxZ if groups come before zones
  if(nest == NEST_GDZ || nest ==  NEST_DGZ || nest == NEST_GZD){
    grid_data->sigt = new SubTVec(NEST_DGZ,
      num_groups, 1, grid_data->num_zones);
  }
  // otherwise, 1xZxG
  else{
    grid_data->sigt = new SubTVec(NEST_DZG,
      num_groups, 1, grid_data->num_zones);
  }
}


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

  error_exit(1);
  return NULL;
}

