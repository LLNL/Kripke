//
// Copyright (c) 2014-19, Lawrence Livermore National Security, LLC
// and Kripke project contributors. See the COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
//

#include <Kripke.h>
#include <Kripke/Arch/LPlusTimes.h>
#include <Kripke/Kernel.h>
#include <Kripke/Timing.h>
#include <Kripke/VarTypes.h>

using namespace Kripke;
using namespace Kripke::Core;

struct LPlusTimesSdom {

  static const std::string KernelName;

  template<typename AL>
  void operator()(AL al, 
                  Kripke::SdomId sdom_id,
                  Set const       &set_dir,
                  Set const       &set_group,
                  Set const       &set_zone,
                  Set const       &set_moment,
                  Field_Moments   &field_phi_out,
                  Field_Flux      &field_rhs,
                  Field_EllPlus   &field_ell_plus) const
  {
    
    using ExecPolicy = typename Kripke::Arch::Policy_LPlusTimes<AL>::ExecPolicy;

    auto sdom_al = getSdomAL(al, sdom_id);

    // Get dimensioning
    int num_directions = set_dir.size(sdom_id);
    int num_groups =     set_group.size(sdom_id);
    int num_moments =    set_moment.size(sdom_id);
    int num_zones =      set_zone.size(sdom_id);

    // Get views
    auto phi_out  = sdom_al.getView(field_phi_out);
    auto rhs      = sdom_al.getView(field_rhs);
    auto ell_plus = sdom_al.getView(field_ell_plus); 

    // Compute:  rhs =  ell_plus * phi_out
    RAJA::kernel<ExecPolicy>(
        camp::make_tuple(
            RAJA::TypedRangeSegment<Direction>(0, num_directions),
            RAJA::TypedRangeSegment<Moment>(0, num_moments),
            RAJA::TypedRangeSegment<Group>(0, num_groups),
            RAJA::TypedRangeSegment<Zone>(0, num_zones) ),
        KRIPKE_LAMBDA (Direction d, Moment nm, Group g, Zone z) {

            rhs(d,g,z) += ell_plus(d, nm) * phi_out(nm, g, z);

        }
    );
  }

};

const std::string LPlusTimesSdom::KernelName =  "LPlusTimes";


void Kripke::Kernel::LPlusTimes(Kripke::Core::DataStore &data_store)
{
  KRIPKE_TIMER(data_store, LPlusTimes);

  Set const &set_dir    = data_store.getVariable<Set>("Set/Direction");
  Set const &set_group  = data_store.getVariable<Set>("Set/Group");
  Set const &set_zone   = data_store.getVariable<Set>("Set/Zone");
  Set const &set_moment = data_store.getVariable<Set>("Set/Moment");

  auto &field_phi_out =   data_store.getVariable<Field_Moments>("phi_out");
  auto &field_rhs =       data_store.getVariable<Field_Flux>("rhs");
  auto &field_ell_plus =  data_store.getVariable<Field_EllPlus>("ell_plus");

  ArchLayoutV al_v = data_store.getVariable<ArchLayout>("al").al_v;

  // Loop over Subdomains
  for (Kripke::SdomId sdom_id : field_rhs.getWorkList()){

    Kripke::dispatch(al_v, LPlusTimesSdom{}, sdom_id,
                     set_dir, set_group, set_zone, set_moment,
                     field_phi_out, field_rhs, field_ell_plus);

  }

}
