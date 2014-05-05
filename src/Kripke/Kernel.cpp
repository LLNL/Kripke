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
      gdset.allocate(grid_data, nest);
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

