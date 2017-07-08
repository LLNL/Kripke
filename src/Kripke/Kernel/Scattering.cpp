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

#include <Kripke/Kernel.h>
#include <Kripke/Grid.h>
#include <Kripke/PartitionSpace.h>
#include <Kripke/SubTVec.h>
#include <Kripke/Timing.h>
#include <Kripke/VarTypes.h>

/**
  Compute scattering source term phi_out from flux moments in phi.
  phi_out(gp,z,nm) = sum_g { sigs(g, n, gp) * phi(g,z,nm) }
*/

void Kripke::Kernel::scattering(Kripke::DataStore &data_store)
{
  KRIPKE_TIMER(data_store, Scattering);

  auto &pspace = data_store.getVariable<Kripke::PartitionSpace>("pspace");

  auto &set_group = data_store.getVariable<Kripke::Set>("Set/Group");

  auto &field_phi =     data_store.getVariable<Kripke::Field_Moments>("phi");
  auto &field_phi_out = data_store.getVariable<Kripke::Field_Moments>("phi_out");
  auto &field_sigs =    data_store.getVariable<Field_SigmaS>("data/sigs");

  Grid_Data *grid_data = &data_store.getVariable<Grid_Data>("grid_data");


  // Zero out source term
  for(auto sdom_id : field_phi_out.getWorkList()){
    auto phi_out_ptr = field_phi_out.getData(sdom_id);
    size_t phi_out_size = field_phi_out.size(sdom_id);
    for(size_t i = 0;i < phi_out_size;++ i){
      phi_out_ptr[i] = 0;
    }
  }

  // Loop over subdomains and compute scattering source
  for(auto sdom_src : field_phi.getWorkList()){
    for(auto sdom_dst : field_phi_out.getWorkList()){

      // Only work on subdomains where src and dst are on the same R subdomain
      size_t r_src = pspace.subdomainToSpace(SPACE_R, sdom_src);
      size_t r_dst = pspace.subdomainToSpace(SPACE_R, sdom_dst);
      if(r_src != r_dst){
        continue;
      }

      // Get glower for src and dst ranges (to index into sigma_s)
      int glower_src = set_group.lower(sdom_src);
      int glower_dst = set_group.lower(sdom_dst);


      // get material mix information
      Subdomain &sdom = grid_data->subdomains[*sdom_src];
      int    const * KRESTRICT zones_to_mixed = &sdom.zones_to_mixed[0];
      int    const * KRESTRICT num_mixed = &sdom.num_mixed[0];
      int    const * KRESTRICT mixed_material = &sdom.mixed_material[0];
      double const * KRESTRICT mixed_fraction = &sdom.mixed_fraction[0];
      int    const * KRESTRICT moment_to_coeff = &grid_data->moment_to_coeff[0];


      auto phi = field_phi.getView(sdom_src);
      auto phi_out = field_phi_out.getView(sdom_dst);
      auto sigs = field_sigs.getView(sdom_src);


      // grab dimensions
      int num_zones = sdom.num_zones;
      int num_groups_src = set_group.size(sdom_src);
      int num_groups_dst = set_group.size(sdom_dst);
      int num_moments = grid_data->total_num_moments;

      for(Moment nm{0};nm < num_moments;++ nm){
        for(Group g{0};g < num_groups_dst;++ g){ // dst group
          for(Group gp{0};gp < num_groups_src;++ gp){ // src group
            for(Zone z{0};z < num_zones;++ z){

              // map nm to n
              Legendre n(moment_to_coeff[*nm]);

              GlobalGroup global_g{*g+glower_dst};
              GlobalGroup global_gp{*gp+glower_src};

              int mix_start = zones_to_mixed[*z];
              int mix_stop = mix_start + num_mixed[*z];

              for(int mix = mix_start;mix < mix_stop;++ mix){
                Material mat{mixed_material[mix]};
                double fraction = mixed_fraction[mix];

                phi_out(nm, g, z) +=
                    sigs(mat, n, global_g, global_gp)
                    * phi(nm, gp, z)
                    * fraction;
              }
            }
          }
        }
      }
    }
  }
}



