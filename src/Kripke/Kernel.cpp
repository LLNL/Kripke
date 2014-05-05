#include<Kripke/Kernel.h>
#include<Kripke/comm.h>
#include<Kripke/grid.h>

#include<Kripke/Kernel/Kernel_3d_GDZ.h>

Kernel::Kernel(){

}
Kernel::~Kernel(){

}

void Kernel::allocateStorage(Grid_Data *grid_data){
  Nesting_Order nest = nesting();
  for(int gs = 0;gs < grid_data->gd_sets.size();++ gs){
    for(int ds = 0;ds < grid_data->gd_sets[gs].size();++ ds){
      Group_Dir_Set &gdset = grid_data->gd_sets[gs][ds];

      delete gdset.psi;
      gdset.psi = new SubTVec(nest,
          gdset.num_groups, gdset.num_directions, grid_data->num_zones);

      delete gdset.rhs;
      gdset.rhs = new SubTVec(nest,
          gdset.num_groups, gdset.num_directions, grid_data->num_zones);

      // allocate sigt  1xGxZ if groups come before zones
      delete gdset.sigt;
      if(nest == NEST_GDZ || nest ==  NEST_DGZ || nest == NEST_GZD){
        gdset.sigt = new SubTVec(NEST_DGZ,
          gdset.num_groups, 1, grid_data->num_zones);
      }
      // otherwise, 1xZxG
      else{
        gdset.sigt = new SubTVec(NEST_DZG,
          gdset.num_groups, 1, grid_data->num_zones);
      }

    }
  }
}


// Factory to create correct kernel object
Kernel *createKernel(Nesting_Order nest, int num_dims){
  if(num_dims == 3){
    switch(nest){
    case NEST_GDZ:
      return new Kernel_3d_GDZ();
    }
  }

  error_exit(1);
  return NULL;
}

