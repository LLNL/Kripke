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

#include<Kripke/Kernel/Kernel_3d_DGZ.h>
#include<Kripke/Grid.h>
#include<Kripke/SubTVec.h>

Nesting_Order Kernel_3d_DGZ::nestingPsi(void) const {
  return NEST_DGZ;
}

Nesting_Order Kernel_3d_DGZ::nestingPhi(void) const {
  return NEST_DGZ;
}

Nesting_Order Kernel_3d_DGZ::nestingSigt(void) const {
  return NEST_DGZ;
}

Nesting_Order Kernel_3d_DGZ::nestingEll(void) const {
  return NEST_ZGD;
}

Nesting_Order Kernel_3d_DGZ::nestingEllPlus(void) const {
  return NEST_ZDG;
}

Nesting_Order Kernel_3d_DGZ::nestingSigs(void) const {
  return NEST_DGZ;
}


void Kernel_3d_DGZ::LTimes(Grid_Data *grid_data) {
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
        for(int gz = 0;gz < num_locgz; ++ gz){
          phi[nm*num_gz + group0*num_zones + gz] += 
            ell[nm*num_local_directions + d] * psi[d*num_locgz + gz];
        }
      }     
    }
  } 
}

void Kernel_3d_DGZ::LPlusTimes(Grid_Data *grid_data) {
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
    int num_groups_zones = num_local_groups*num_zones;
    
    // Get pointers
    double const * KRESTRICT phi_out = sdom.phi_out->ptr() + group0*num_zones;
    double const * KRESTRICT ell_plus = sdom.ell_plus->ptr();
    double       * KRESTRICT rhs = sdom.rhs->ptr();

    for (int d = 0; d < num_local_directions; d++) {      
      for(int nm = 0;nm < num_moments;++nm){
#ifdef KRIPKE_USE_OPENMP
#pragma omp parallel for
#endif
        for(int gz = 0;gz < num_groups_zones; ++ gz){
          rhs[d*num_groups_zones + gz] +=
            ell_plus[d*num_moments + nm] * phi_out[nm*num_groups*num_zones + gz];
        }
      }
    }
  }
}

/**
  Compute scattering source term phi_out from flux moments in phi.
  phi_out(gp,z,nm) = sum_g { sigs(g, n, gp) * phi(g,z,nm) }
*/
void Kernel_3d_DGZ::scattering(Grid_Data *grid_data){
  // Zero zone sets
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
      for(int g = 0;g < num_groups;++ g){      
        for(int gp = 0;gp < num_groups;++ gp){
#ifdef KRIPKE_USE_OPENMP
#pragma omp parallel for
#endif
          for(int zone = 0;zone < num_zones;++ zone){            
            int mix_start = zones_to_mixed[zone];
            int mix_stop = mix_start + num_mixed[zone];

            for(int mix = mix_start;mix < mix_stop;++ mix){
              int material = mixed_material[mix];
              double fraction = mixed_fraction[mix];              

              phi_out[nm*num_gz + gp*num_zones + zone] += 
                sigs[n*3*num_groups*num_groups + g*3*num_groups + gp*3 + material] * 
                phi[nm*num_gz + g*num_zones + zone] * fraction;
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
void Kernel_3d_DGZ::source(Grid_Data *grid_data){
  // Loop over zoneset subdomains
  for(int zs = 0;zs < grid_data->num_zone_sets;++ zs){
    // get the phi and phi out references
    SubTVec &phi_out = *grid_data->phi_out[zs];

    // get material mix information
    int sdom_id = grid_data->zs_to_sdomid[zs];
    Subdomain &sdom = grid_data->subdomains[sdom_id];
    int    const * KRESTRICT mixed_to_zones = &sdom.mixed_to_zones[0];
    int    const * KRESTRICT mixed_material = &sdom.mixed_material[0];
    double const * KRESTRICT mixed_fraction = &sdom.mixed_fraction[0];
    double       * KRESTRICT phi_out_nm0 = phi_out.ptr();

    // grab dimensions
    int num_mixed = sdom.mixed_to_zones.size();
    int num_zones = sdom.num_zones;
    int num_groups = phi_out.groups;

#ifdef KRIPKE_USE_OPENMP
#pragma omp parallel for
#endif
    for(int g = 0;g < num_groups;++ g){
      double       * KRESTRICT phi_out_nm0_g = phi_out_nm0 + g*num_zones;
      
      for(int mix = 0;mix < num_mixed;++ mix){
        int zone = mixed_to_zones[mix];
        int material = mixed_material[mix];
        double fraction = mixed_fraction[mix];

        if(material == 0){
          phi_out_nm0_g[zone] += 1.0 * fraction;
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

void Kernel_3d_DGZ::sweep(Subdomain *sdom) {
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
  
  int num_gz = num_groups * num_zones;
  int num_gz_i = local_jmax * local_kmax * num_groups;
  int num_gz_j = local_imax * local_kmax * num_groups;
  int num_gz_k = local_imax * local_jmax * num_groups;
  int num_z_i = local_jmax * local_kmax;
  int num_z_j = local_imax * local_kmax;
  int num_z_k = local_imax * local_jmax;

  // All directions have same id,jd,kd, since these are all one Direction Set
  // So pull that information out now
  Grid_Sweep_Block const &extent = sdom->sweep_block;

#ifdef KRIPKE_USE_OPENMP
#pragma omp parallel for
#endif
  for (int d = 0; d < num_directions; ++d) {
    for (int g = 0; g < num_groups; ++g) {
      for (int k = extent.start_k; k != extent.end_k; k += extent.inc_k) {               
        for (int j = extent.start_j; j != extent.end_j; j += extent.inc_j) {          
          for (int i = extent.start_i; i != extent.end_i; i += extent.inc_i) {
            double const xcos_dxi = 2.0 * direction[d].xcos / dx[i + 1];
            double const ycos_dyj = 2.0 * direction[d].ycos / dy[j + 1];
            double const zcos_dzk = 2.0 * direction[d].zcos / dz[k + 1];
            
            int const z_idx = Zonal_INDEX(i, j, k);            
            int const lf_idx = d*num_gz_i + g*num_z_i + I_PLANE_INDEX(j, k);
            int const fr_idx = d*num_gz_j + g*num_z_j + J_PLANE_INDEX(i, k);
            int const bo_idx = d*num_gz_k + g*num_z_k + K_PLANE_INDEX(i, j);

            /* Calculate new zonal flux */
            double const psi_d_g_z = (
                  rhs[d*num_gz + g*num_zones + z_idx]
                + psi_lf[lf_idx] * xcos_dxi
                + psi_fr[fr_idx] * ycos_dyj
                + psi_bo[bo_idx] * zcos_dzk)
                / (xcos_dxi + ycos_dyj + zcos_dzk + sigt[g*num_zones + z_idx]);

            psi[d*num_gz + g*num_zones + z_idx] = psi_d_g_z;
            
            /* Apply diamond-difference relationships */
            psi_lf[lf_idx] = 2.0 * psi_d_g_z - psi_lf[lf_idx];
            psi_fr[fr_idx] = 2.0 * psi_d_g_z - psi_fr[fr_idx];
            psi_bo[bo_idx] = 2.0 * psi_d_g_z - psi_bo[bo_idx];
          }
        }
      }
    } // group
  } // direction

}


