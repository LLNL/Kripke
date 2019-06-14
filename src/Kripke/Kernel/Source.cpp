//
// Copyright (c) 2014-19, Lawrence Livermore National Security, LLC
// and Kripke project contributors. See the COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
//

#include <Kripke/Kernel.h>

#include <Kripke.h>
#include <Kripke/Arch/Source.h>
#include <Kripke/Timing.h>
#include <Kripke/VarTypes.h>

using namespace Kripke;
using namespace Kripke::Core;

/**
 * Add an isotropic source, with flux of 1, to every zone with Region 1
 * (or material 0).
 *
 * Since it's isotropic, we're just adding this to nm=0.
 */
struct SourceSdom {

  static const std::string KernelName;

  template<typename AL>
  RAJA_INLINE
  void operator()(AL al, 
                  Kripke::SdomId          sdom_id,
                  Set const               &set_group,
                  Set const               &set_mixelem,
                  Field_Moments           &field_phi_out,
                  Field_MixElem2Zone      &field_mixed_to_zone,
                  Field_MixElem2Material  &field_mixed_to_material,
                  Field_MixElem2Double    &field_mixed_to_fraction,
                  double                  source_strength) const
  {

    using ExecPolicy = typename Kripke::Arch::Policy_Source<AL>::ExecPolicy;

    auto sdom_al = getSdomAL(al, sdom_id);

    // Source term is isotropic
    Moment nm{0};

    auto phi_out = sdom_al.getView(field_phi_out);

    auto mixelem_to_zone     = sdom_al.getView(field_mixed_to_zone);
    auto mixelem_to_material = sdom_al.getView(field_mixed_to_material);
    auto mixelem_to_fraction = sdom_al.getView(field_mixed_to_fraction);

    int num_mixed  = set_mixelem.size(sdom_id);
    int num_groups = set_group.size(sdom_id);


    // Compute:  phi =  ell * psi
    RAJA::kernel<ExecPolicy>(
        camp::make_tuple(
            RAJA::TypedRangeSegment<Group>(0, num_groups),
            RAJA::TypedRangeSegment<MixElem>(0, num_mixed) ),
        KRIPKE_LAMBDA (Group g, MixElem mix) {

            Material material = mixelem_to_material(mix);

            if(material == 0){
              Zone z = mixelem_to_zone(mix);
              double fraction = mixelem_to_fraction(mix);

              phi_out(nm, g, z) += source_strength * fraction;
            }

        }
    );

  }
};

const std::string SourceSdom::KernelName = "Source";


void Kripke::Kernel::source(DataStore &data_store)
{
  KRIPKE_TIMER(data_store, Source);

  ArchLayoutV al_v = data_store.getVariable<ArchLayout>("al").al_v;

  auto &set_group   = data_store.getVariable<Set>("Set/Group");
  auto &set_mixelem = data_store.getVariable<Set>("Set/MixElem");

  auto &field_phi_out = data_store.getVariable<Kripke::Field_Moments>("phi_out");

  auto &field_mixed_to_zone     = data_store.getVariable<Field_MixElem2Zone>("mixelem_to_zone");
  auto &field_mixed_to_material = data_store.getVariable<Field_MixElem2Material>("mixelem_to_material");
  auto &field_mixed_to_fraction = data_store.getVariable<Field_MixElem2Double>("mixelem_to_fraction");

  double source_strength = 1.0;


  // Loop over zoneset subdomains
  for(auto sdom_id : field_phi_out.getWorkList()){

    Kripke::dispatch(al_v, SourceSdom{}, sdom_id,
                     set_group, set_mixelem,
                     field_phi_out,
                     field_mixed_to_zone,
                     field_mixed_to_material,
                     field_mixed_to_fraction,
                     source_strength);

  }


}
