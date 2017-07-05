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
#include<Kripke/Timing.h>

/**
  Compute scattering source term phi_out from flux moments in phi.
  phi_out(gp,z,nm) = sum_g { sigs(g, n, gp) * phi(g,z,nm) }
*/

void Kripke::Kernel::scattering(Kripke::DataStore &data_store)
{
  KRIPKE_TIMER(data_store, Scattering);

  Grid_Data *grid_data = &data_store.getVariable<Grid_Data>("grid_data");

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

    // Zero out source terms
    grid_data->phi_out[zs]->clear(0.0);

    // grab dimensions
    //int num_mixed = sdom.mixed_to_zones.size();
    int num_zones = sdom.num_zones;
    int num_groups = grid_data->phi_out[zs]->groups;
    int num_moments = grid_data->total_num_moments;
    int num_gz = num_groups*num_zones;

    for(int nm = 0;nm < num_moments;++ nm){
      // map nm to n
      int n = moment_to_coeff[nm];
      double const * KRESTRICT sigs_n = sigs + n*3*num_groups*num_groups;
      double const * KRESTRICT phi_nm = phi + nm*num_gz;
      double       * KRESTRICT phi_out_nm = phi_out + nm*num_gz;

      for(int g = 0;g < num_groups;++ g){      
        double const * KRESTRICT sigs_n_g = sigs_n + g*3*num_groups;
        double const * KRESTRICT phi_nm_g = phi_nm + g*num_zones;
                
        for(int gp = 0;gp < num_groups;++ gp){
          double const * KRESTRICT sigs_n_g_gp = sigs_n_g + gp*3;
          double       * KRESTRICT phi_out_nm_gp = phi_out_nm + gp*num_zones;

          for(int zone = 0;zone < num_zones;++ zone){
            double phi_out_nm_gp_z = 0.0;
            int mix_start = zones_to_mixed[zone];
            int mix_stop = mix_start + num_mixed[zone];

            for(int mix = mix_start;mix < mix_stop;++ mix){
              int material = mixed_material[mix];
              double fraction = mixed_fraction[mix];

              phi_out_nm_gp_z += sigs_n_g_gp[material] * phi_nm_g[zone] * fraction;
            }
            phi_out_nm_gp[zone] += phi_out_nm_gp_z;
          }
        }        
      }
    }
  }
}



