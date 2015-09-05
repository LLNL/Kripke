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

#include<Kripke/Kernel/Kernel_3d_DZG.h>
#include<Kripke/Grid.h>
#include<Kripke/SubTVec.h>

Nesting_Order Kernel_3d_DZG::nestingPsi(void) const {
  return NEST_DZG;
}

Nesting_Order Kernel_3d_DZG::nestingPhi(void) const {
  return NEST_DZG;
}

Nesting_Order Kernel_3d_DZG::nestingSigt(void) const {
  return NEST_DZG;
}

Nesting_Order Kernel_3d_DZG::nestingEll(void) const {
  return NEST_ZGD;
}

Nesting_Order Kernel_3d_DZG::nestingEllPlus(void) const {
  return NEST_ZDG;
}

Nesting_Order Kernel_3d_DZG::nestingSigs(void) const {
  return NEST_DZG;
}


void Kernel_3d_DZG::LTimes(Grid_Data *grid_data) {
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
    int num_gz = num_groups*num_zones;
    int num_locgz = num_local_groups*num_zones;
    
    // Get pointers
    double const * KRESTRICT ell = sdom.ell->ptr();
    double const * KRESTRICT psi = sdom.psi->ptr();
    double       * KRESTRICT phi = sdom.phi->ptr();
    
    for(int nm = 0;nm < num_moments;++nm){     
      for (int d = 0; d < num_local_directions; d++) {
#ifdef KRIPKE_USE_OPENMP
#pragma omp parallel for
#endif
        for (int z = 0;z < num_zones;++ z){                         
          for(int g = 0;g < num_local_groups; ++ g){  
            phi[nm*num_gz + z*num_groups + group0 + g] += 
              ell[nm*num_local_directions + d] * 
              psi[d*num_locgz + z*num_local_groups + g];
          }
        }
      }
    }
  } 
}

void Kernel_3d_DZG::LPlusTimes(Grid_Data *grid_data) {
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
    int num_groups = sdom.phi_out->groups;
    int num_local_groups = sdom.num_groups;
    int group0 = sdom.group0;
    int num_local_directions = sdom.num_directions;
    int num_gz = num_groups*num_zones;
    int num_locgz = num_local_groups*num_zones;

    // Get pointers
    double const * KRESTRICT phi_out = sdom.phi_out->ptr() + group0;
    double const * KRESTRICT ell_plus = sdom.ell_plus->ptr();
    double       * KRESTRICT rhs = sdom.rhs->ptr();

    for (int d = 0; d < num_local_directions; d++) {
      for(int nm = 0;nm < num_moments;++nm){
#ifdef KRIPKE_USE_OPENMP
#pragma omp parallel for
#endif
        for(int z = 0;z < num_zones;++ z){
          for(int g = 0;g < num_local_groups;++ g){
            rhs[d*num_locgz + z*num_local_groups + g] += 
              ell_plus[d*num_moments + nm] * 
              phi_out[nm*num_gz + z*num_groups + g];
          }
        }
      }
    }
  }
}


/**
  Compute scattering source term phi_out from flux moments in phi.
  phi_out(gp,z,nm) = sum_g { sigs(g, n, gp) * phi(g,z,nm) }

*/
void Kernel_3d_DZG::scattering(Grid_Data *grid_data){
  // Zero out source terms
  for(int zs = 0;zs < grid_data->num_zone_sets;++ zs){
    grid_data->phi_out[zs]->clear(0.0);
  }
  
  // Loop over zoneset subdomains
  for(int zs = 0;zs < grid_data->num_zone_sets;++ zs){
    // get material mix information
    int sdom_id = grid_data->zs_to_sdomid[zs];
    Subdomain &sdom = grid_data->subdomains[sdom_id];
    int    const * KRESTRICT zones_to_mixed = &sdom.zones_to_mixed[0];
    int    const * KRESTRICT num_mixed = &sdom.num_mixed[0];
    int    const * KRESTRICT mixed_material = &sdom.mixed_material[0];
    double const * KRESTRICT mixed_fraction = &sdom.mixed_fraction[0];
    double const * KRESTRICT sigs = grid_data->sigs->ptr(); 

    int    const * KRESTRICT moment_to_coeff = &grid_data->moment_to_coeff[0];
    double const * KRESTRICT phi = grid_data->phi[zs]->ptr();
    double       * KRESTRICT phi_out = grid_data->phi_out[zs]->ptr();

    // grab dimensions
    int num_zones = sdom.num_zones;
    int num_groups = grid_data->phi_out[zs]->groups;
    int num_moments = grid_data->total_num_moments;
    int num_gz = num_groups*num_zones;

    for(int nm = 0;nm < num_moments;++ nm){
      int n = moment_to_coeff[nm];
#ifdef KRIPKE_USE_OPENMP
#pragma omp parallel for
#endif
      for(int zone = 0;zone < num_zones;++ zone){
        int mix_start = zones_to_mixed[zone];
        int mix_stop = mix_start + num_mixed[zone];

        for(int mix = mix_start;mix < mix_stop;++ mix){
          int material = mixed_material[mix];
          double fraction = mixed_fraction[mix];       
          for(int g = 0;g < num_groups;++ g){                          
            for(int gp = 0;gp < num_groups;++ gp){           
                                                      
              phi_out[nm*num_gz + zone*num_groups + gp] += 
                sigs[n*3*num_groups*num_groups + material*num_groups*num_groups + g*num_groups + gp] * 
                phi[nm*num_gz + zone*num_groups + g] * 
                fraction;
            }
          }        
        }
      }
    }
  }
}

/**
 * Add an isotropic source, with flux of 1, to every zone with Region 1
 * (or material 0).
 *
 * Since it's isotropic, we're just adding this to nm=0.
 */
void Kernel_3d_DZG::source(Grid_Data *grid_data){
  // Loop over zoneset subdomains
  for(int zs = 0;zs < grid_data->num_zone_sets;++ zs){
    // get material mix information
    int sdom_id = grid_data->zs_to_sdomid[zs];
    Subdomain &sdom = grid_data->subdomains[sdom_id];
    int    const * KRESTRICT zones_to_mixed = &sdom.zones_to_mixed[0];
    int    const * KRESTRICT num_mixed = &sdom.num_mixed[0];
    int    const * KRESTRICT mixed_material = &sdom.mixed_material[0];
    double const * KRESTRICT mixed_fraction = &sdom.mixed_fraction[0];
    double       * KRESTRICT phi_out = grid_data->phi_out[zs]->ptr();

    // grab dimensions
    int num_zones = sdom.num_zones;
    int num_groups = grid_data->phi_out[zs]->groups;
    int num_moments = grid_data->total_num_moments;
    
#ifdef KRIPKE_USE_OPENMP
#pragma omp parallel for
#endif
    for(int zone = 0;zone < num_zones;++ zone){
      int mix_start = zones_to_mixed[zone];
      int mix_stop = mix_start + num_mixed[zone];

      for(int mix = mix_start;mix < mix_stop;++ mix){
        int material = mixed_material[mix];
        double fraction = mixed_fraction[mix];
        
        for(int g = 0;g < num_groups;++ g){
          if(material == 0){
          
            phi_out[zone*num_groups + g] += 1.0 * fraction;
          }
        }
      }
    }
  }
}



// Macros for offsets with fluxes on cell faces 
#define I_PLANE_INDEX(j, k) ((k)*(local_jmax) + (j))
#define J_PLANE_INDEX(i, k) ((k)*(local_imax) + (i))
#define K_PLANE_INDEX(i, j) ((j)*(local_imax) + (i))
#define Zonal_INDEX(i, j, k) ((i) + (local_imax)*(j) \
  + (local_imax)*(local_jmax)*(k))
  
void Kernel_3d_DZG::sweep(Subdomain *sdom) {
  int num_directions = sdom->num_directions;
  int num_groups = sdom->num_groups;
  int num_zones = sdom->num_zones;

  Directions *direction = sdom->directions;

  int local_imax = sdom->nzones[0];
  int local_jmax = sdom->nzones[1];
  int local_kmax = sdom->nzones[2];

  double const * KRESTRICT dx = &sdom->deltas[0][0];
  double const * KRESTRICT dy = &sdom->deltas[1][0];
  double const * KRESTRICT dz = &sdom->deltas[2][0];
  
  double const * KRESTRICT sigt = sdom->sigt->ptr();
  double       * KRESTRICT psi  = sdom->psi->ptr();   
  double const * KRESTRICT rhs  = sdom->rhs->ptr();

  double * KRESTRICT psi_lf = sdom->plane_data[0]->ptr();
  double * KRESTRICT psi_fr = sdom->plane_data[1]->ptr();
  double * KRESTRICT psi_bo = sdom->plane_data[2]->ptr();
  
  int num_zg = num_zones * num_groups;
  int num_zg_i = local_jmax * local_kmax * num_groups;
  int num_zg_j = local_imax * local_kmax * num_groups;
  int num_zg_k = local_imax * local_jmax * num_groups;

  // All directions have same id,jd,kd, since these are all one Direction Set
  // So pull that information out now
  Grid_Sweep_Block const &extent = sdom->sweep_block;

#ifdef KRIPKE_USE_OPENMP
#pragma omp parallel for
#endif
  for (int d = 0; d < num_directions; ++d) {
    for (int k = extent.start_k; k != extent.end_k; k += extent.inc_k) {            
      for (int j = extent.start_j; j != extent.end_j; j += extent.inc_j) {                
        for (int i = extent.start_i; i != extent.end_i; i += extent.inc_i) {
          for (int g = 0; g < num_groups; ++g) {                      
            double const xcos_dxi = 2.0 * direction[d].xcos / dx[i + 1];
            double const ycos_dyj = 2.0 * direction[d].ycos / dy[j + 1];
            double const zcos_dzk = 2.0 * direction[d].zcos / dz[k + 1];
            
            int const z = Zonal_INDEX(i, j, k);
            int const lf_idx = d*num_zg_i + I_PLANE_INDEX(j, k)*num_groups + g;
            int const fr_idx = d*num_zg_j + J_PLANE_INDEX(i, k)*num_groups + g;
            int const bo_idx = d*num_zg_k + K_PLANE_INDEX(i, j)*num_groups + g;
            
            // Calculate new zonal flux 
            double const psi_d_z_g = (
                rhs[d*num_zg + z*num_groups + g]
                + psi_lf[lf_idx] * xcos_dxi
                + psi_fr[fr_idx] * ycos_dyj
                + psi_bo[bo_idx] * zcos_dzk)
                / (xcos_dxi + ycos_dyj + zcos_dzk + sigt[z*num_groups + g]);

            psi[d*num_zg + z*num_groups + g] = psi_d_z_g;

            // Apply diamond-difference relationships 
            psi_lf[lf_idx] = 2.0 * psi_d_z_g - psi_lf[lf_idx];
            psi_fr[fr_idx] = 2.0 * psi_d_z_g - psi_fr[fr_idx];
            psi_bo[bo_idx] = 2.0 * psi_d_z_g - psi_bo[bo_idx];
          }
        }
      }
    }
  }
}

