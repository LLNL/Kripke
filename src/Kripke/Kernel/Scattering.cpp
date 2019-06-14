//
// Copyright (c) 2014-19, Lawrence Livermore National Security, LLC
// and Kripke project contributors. See the COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
//

#include <Kripke/Kernel.h>

#include <Kripke.h>
#include <Kripke/Arch/Scattering.h>
#include <Kripke/Core/PartitionSpace.h>
#include <Kripke/Timing.h>
#include <Kripke/VarTypes.h>

using namespace Kripke;
using namespace Kripke::Core;

/**
  Compute scattering source term phi_out from flux moments in phi.
  phi_out(gp,z,nm) = sum_g { sigs(g, n, gp) * phi(g,z,nm) }
*/

struct ScatteringSdom {

  static const std::string KernelName;

  template<typename AL>
  RAJA_INLINE
  void operator()(AL al, 
                  Kripke::SdomId          sdom_src,
                  Kripke::SdomId          sdom_dst,
                  Set const               &set_group,
                  Set const               &set_zone,
                  Set const               &set_moment,
                  Field_Moments           &field_phi,
                  Field_Moments           &field_phi_out,
                  Field_SigmaS            &field_sigs,
                  Field_Zone2MixElem      &field_zone_to_mixelem,
                  Field_Zone2Int          &field_zone_to_num_mixelem,
                  Field_MixElem2Material  &field_mixelem_to_material,
                  Field_MixElem2Double    &field_mixelem_to_fraction,
                  Field_Moment2Legendre   &field_moment_to_legendre) const
  {

    using ExecPolicy = typename Kripke::Arch::Policy_Scattering<AL>::ExecPolicy;

    auto sdom_al = getSdomAL(al, sdom_src);

    // Get glower for src and dst ranges (to index into sigma_s)
    int glower_src = set_group.lower(sdom_src);
    int glower_dst = set_group.lower(sdom_dst);


    // get material mix information
    auto moment_to_legendre = sdom_al.getView(field_moment_to_legendre);

    auto phi     = sdom_al.getView(field_phi);
    auto phi_out = sdom_al.getView(field_phi_out, sdom_dst);
    auto sigs    = sdom_al.getView(field_sigs);
    
    auto zone_to_mixelem     = sdom_al.getView(field_zone_to_mixelem);
    auto zone_to_num_mixelem = sdom_al.getView(field_zone_to_num_mixelem);
    auto mixelem_to_material = sdom_al.getView(field_mixelem_to_material);
    auto mixelem_to_fraction = sdom_al.getView(field_mixelem_to_fraction);
    
    // grab dimensions
    int num_zones =      set_zone.size(sdom_src);
    int num_groups_src = set_group.size(sdom_src);
    int num_groups_dst = set_group.size(sdom_dst);
    int num_moments =    set_moment.size(sdom_dst);

    RAJA::kernel<ExecPolicy>(
        camp::make_tuple(
            RAJA::TypedRangeSegment<Moment>(0, num_moments),
            RAJA::TypedRangeSegment<Group>(0, num_groups_dst),
            RAJA::TypedRangeSegment<Group>(0, num_groups_src),
            RAJA::TypedRangeSegment<Zone>(0, num_zones) ),
        KRIPKE_LAMBDA (Moment nm, Group g, Group gp, Zone z) {

            // map nm to n
            Legendre n = moment_to_legendre(nm);

            GlobalGroup global_g{*g+glower_dst};
            GlobalGroup global_gp{*gp+glower_src};

            MixElem mix_start = zone_to_mixelem(z);
            MixElem mix_stop = mix_start + zone_to_num_mixelem(z);

            double sigs_z = 0.0;
            for(MixElem mix = mix_start;mix < mix_stop;++ mix){
              Material mat = mixelem_to_material(mix);
              double fraction = mixelem_to_fraction(mix);

              sigs_z += sigs(mat, n, global_g, global_gp) * fraction;
            }
            phi_out(nm, g, z) += sigs_z * phi(nm, gp, z);
        }
    );
  }

};

const std::string ScatteringSdom::KernelName = "Scattering";


/**
  Compute scattering source term phi_out from flux moments in phi.
  phi_out(gp,z,nm) = sum_g { sigs(g, n, gp) * phi(g,z,nm) }
*/

void Kripke::Kernel::scattering(Kripke::Core::DataStore &data_store)
{
  KRIPKE_TIMER(data_store, Scattering);

  ArchLayoutV al_v = data_store.getVariable<ArchLayout>("al").al_v;

  auto &pspace = data_store.getVariable<Kripke::Core::PartitionSpace>("pspace");

  auto &set_group  = data_store.getVariable<Kripke::Core::Set>("Set/Group");
  auto &set_moment = data_store.getVariable<Kripke::Core::Set>("Set/Moment");
  auto &set_zone   = data_store.getVariable<Kripke::Core::Set>("Set/Zone");

  auto &field_phi     = data_store.getVariable<Kripke::Field_Moments>("phi");
  auto &field_phi_out = data_store.getVariable<Kripke::Field_Moments>("phi_out");
  auto &field_sigs    = data_store.getVariable<Field_SigmaS>("data/sigs");

  auto &field_zone_to_mixelem     = data_store.getVariable<Field_Zone2MixElem>("zone_to_mixelem");
  auto &field_zone_to_num_mixelem = data_store.getVariable<Field_Zone2Int>("zone_to_num_mixelem");
  auto &field_mixed_to_material = data_store.getVariable<Field_MixElem2Material>("mixelem_to_material");
  auto &field_mixed_to_fraction = data_store.getVariable<Field_MixElem2Double>("mixelem_to_fraction");

  auto &field_moment_to_legendre = data_store.getVariable<Field_Moment2Legendre>("moment_to_legendre");


  // Loop over subdomains and compute scattering source
  for(auto sdom_src : field_phi.getWorkList()){
    for(auto sdom_dst : field_phi_out.getWorkList()){

      // Only work on subdomains where src and dst are on the same R subdomain
      size_t r_src = pspace.subdomainToSpace(SPACE_R, sdom_src);
      size_t r_dst = pspace.subdomainToSpace(SPACE_R, sdom_dst);
      if(r_src != r_dst){
        continue;
      }

      Kripke::dispatch(al_v, ScatteringSdom{}, sdom_src,
                       sdom_dst,
                       set_group, set_zone, set_moment,
                       field_phi, field_phi_out, field_sigs,
                       field_zone_to_mixelem,
                       field_zone_to_num_mixelem,
                       field_mixed_to_material,
                       field_mixed_to_fraction,
                       field_moment_to_legendre);



    }

  }


}


