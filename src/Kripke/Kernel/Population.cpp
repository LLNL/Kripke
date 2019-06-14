//
// Copyright (c) 2014-19, Lawrence Livermore National Security, LLC
// and Kripke project contributors. See the COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
//

#include <Kripke/Kernel.h>

#include <Kripke.h>
#include <Kripke/Arch/Population.h>
#include <Kripke/Core/Comm.h>
#include <Kripke/Timing.h>
#include <Kripke/VarTypes.h>

using namespace Kripke;
using namespace Kripke::Core;


struct PopulationSdom {

  static const std::string KernelName;

  template<typename AL>
  void operator()(AL al, 
                  Kripke::SdomId sdom_id,
                  Set const               &set_dir,
                  Set const               &set_group,
                  Set const               &set_zone,
                  Field_Flux              &field_psi,
                  Field_Direction2Double  &field_w,
                  Field_Zone2Double       &field_volume,
                  double                  *part_ptr) const
  {
    using Policy = Kripke::Arch::Policy_Population<AL>;
    using ReducePolicy = typename Policy::ReducePolicy;
    using ExecPolicy = typename Policy::ExecPolicy;

    auto sdom_al = getSdomAL(al, sdom_id);

    int num_directions = set_dir.size(sdom_id);
    int num_groups =     set_group.size(sdom_id);
    int num_zones =      set_zone.size(sdom_id);

    auto psi    = sdom_al.getView(field_psi);
    auto w      = sdom_al.getView(field_w);
    auto volume = sdom_al.getView(field_volume);
    
    RAJA::ReduceSum<ReducePolicy, double> part_red(0.0);

    RAJA::kernel<ExecPolicy>(
        camp::make_tuple(
            RAJA::TypedRangeSegment<Direction>(0, num_directions),
            RAJA::TypedRangeSegment<Group>(0, num_groups),
            RAJA::TypedRangeSegment<Zone>(0, num_zones) ),
        KRIPKE_LAMBDA (Direction d, Group g, Zone z) {

          part_red += w(d) * psi(d,g,z) * volume(z);

        }
    );

    *part_ptr += (double)part_red;
  }

};

const std::string PopulationSdom::KernelName = "Population";


/**
 * Returns the integral of Psi over all phase-space, to look at convergence
 */
double Kripke::Kernel::population(Kripke::Core::DataStore &data_store)
{
  KRIPKE_TIMER(data_store, Population);

  ArchLayoutV al_v = data_store.getVariable<ArchLayout>("al").al_v; 

  Set const &set_dir    = data_store.getVariable<Set>("Set/Direction");
  Set const &set_group  = data_store.getVariable<Set>("Set/Group");
  Set const &set_zone   = data_store.getVariable<Set>("Set/Zone");

  auto &field_psi =       data_store.getVariable<Field_Flux>("psi");
  auto &field_w =         data_store.getVariable<Field_Direction2Double>("quadrature/w");
  auto &field_volume =    data_store.getVariable<Field_Zone2Double>("volume");

  // sum up particles for psi and rhs
  double part = 0.0;
  for (Kripke::SdomId sdom_id : field_psi.getWorkList()){

    Kripke::dispatch(al_v, PopulationSdom{}, sdom_id,
                     set_dir, set_group, set_zone,
                     field_psi, field_w, field_volume,
                     &part);
  }

  // reduce
  auto const &comm = data_store.getVariable<Kripke::Core::Comm>("comm");
  return comm.allReduceSumDouble(part);
}

