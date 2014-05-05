#include<Kripke/Kernel/Kernel_GZ_Base.h>
#include<Kripke/user_data.h>
#include<Kripke/SubTVec.h>

Kernel_GZ_Base::Kernel_GZ_Base(){

}

Kernel_GZ_Base::~Kernel_GZ_Base(){

}

// Computational Kernels
void Kernel_GZ_Base::evalSigmaTot(User_Data *user_data, Group_Dir_Set *ga_set){
  int num_groups = ga_set->num_groups;
  int num_zones = user_data->grid_data->num_zones;
  double **sigt = ga_set->sigt->data[0];
  double *grp_sigt = &user_data->sigma_tot[0];

  for(int g = 0;g < num_groups;++ g){
    double *sigt_g = sigt[g];
    for(int z = 0;z < num_zones;++ z){
      sigt_g[z] = grp_sigt[z];
    }
  }
}


