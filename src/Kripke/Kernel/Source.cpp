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
 * Add an isotropic source, with flux of 1, to every zone with Region 1
 * (or material 0).
 *
 * Since it's isotropic, we're just adding this to nm=0.
 */
void Kripke::Kernel::source(Kripke::DataStore &data_store)
{
  KRIPKE_TIMER(data_store, Source);

  Grid_Data *grid_data = &data_store.getVariable<Grid_Data>("grid_data");

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

