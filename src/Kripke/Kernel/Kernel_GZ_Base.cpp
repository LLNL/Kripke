#include<Kripke/Kernel/Kernel_GZ_Base.h>
#include<Kripke/User_Data.h>
#include<Kripke/SubTVec.h>

Kernel_GZ_Base::Kernel_GZ_Base(){

}

Kernel_GZ_Base::~Kernel_GZ_Base(){

}

// Computational Kernels
void Kernel_GZ_Base::evalSigmaTot(User_Data *user_data, Group_Dir_Set *ga_set){
  int num_groups = ga_set->num_groups;
  int num_zones = user_data->grid_data->num_zones;
  double *sigt_table = &user_data->sigma_tot[0];

  for(int g = 0;g < num_groups;++ g){
    double *sigt_g = ga_set->sigt->ptr(g, 0, 0);
    double g_sigt_table = sigt_table[g];

    for(int z = 0;z < num_zones;++ z){
      sigt_g[z] = g_sigt_table;
    }
  }
}


