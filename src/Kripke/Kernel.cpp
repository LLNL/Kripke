#include<Kripke/Kernel.h>
#include<Kripke/comm.h>
#include<Kripke/grid.h>
#include<Kripke/user_data.h>
#include<Kripke/LMat.h>
#include<Kripke/SubTVec.h>

#include<Kripke/Kernel/Kernel_3d_GDZ.h>

Kernel::Kernel(){

}
Kernel::~Kernel(){

}

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

  grid_data->ell = new LMat(NEST_NMD, num_dims, num_moments, total_dirs);
  grid_data->ell_plus = new LMat(NEST_DNM, num_dims, num_moments, total_dirs);

  int total_moments = grid_data->ell->total_moments;
  grid_data->phi = new SubTVec(nestingPhi(), num_groups, total_moments, num_zones);
  grid_data->phi_out = new SubTVec(nestingPhi(), num_groups, total_moments, num_zones);
}


void Kernel::evalSigmaS(Grid_Data *grid_data, int n, int g, int gp){

  int num_zones = grid_data->num_zones;
  double *sig_s = &grid_data->sig_s[0];

  for(int z = 0;z < num_zones;++ z){
    sig_s[z] = 0.0;
  }
}


void Kernel::scattering(Grid_Data *grid_data){
  int num_moments = grid_data->num_moments;
  int num_groups = grid_data->phi->groups;
  int num_zones = grid_data->num_zones;

  double ***phi_in = grid_data->phi->data;
  double ***phi_out = grid_data->phi_out->data;

  // Loop over destination group
  for(int gp=0; gp < num_groups; gp++){

    // Begin loop over scattering moments
    for(int n=0; n < num_moments; n++){

      int num_m = grid_data->ell->numM(n);

      // Loop over source group
      int m0 = 0;
      for(int g=0; g < num_groups; g++){

        // Evaluate sigs  for this (n,g,gp) triplet
        evalSigmaS(grid_data, n, g, gp);

        // Get variables
        double *sig_s = &grid_data->sig_s[0];
        double **phi_in_g = phi_in[g];
        double **phi_out_g = phi_out[g];

        for(int m=0; m < num_m; m++){
          double * __restrict__ phi_out_g_nm = phi_out_g[m+m0];
          double * __restrict__ phi_in_g_nm = phi_in_g[m+m0];

          for(int zone=0; zone<num_zones; zone++){
            phi_out_g_nm[zone] += sig_s[zone]*phi_in_g_nm[zone];
          }

        } // m
      } // g

      m0 += num_m;
    } // n
  } // gp
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

