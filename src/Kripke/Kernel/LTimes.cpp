//
// Copyright (c) 2014-19, Lawrence Livermore National Security, LLC
// and Kripke project contributors. See the COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
//

#include <Kripke.h>
#include <Kripke/Arch/LTimes.h>
#include <Kripke/Core/VarLayout.h>
#include <Kripke/Kernel.h>
#include <Kripke/Timing.h>
#include <Kripke/VarTypes.h>

#include<utility>

using namespace Kripke;
using namespace Kripke::Core;


struct LTimesSdom {

  static const std::string KernelName;

  template<typename AL>
  RAJA_INLINE
  void operator()(AL al, 
                  Kripke::SdomId sdom_id,
                  Set const       &set_dir,
                  Set const       &set_group,
                  Set const       &set_zone,
                  Set const       &set_moment,
                  Field_Flux      &field_psi,
                  Field_Moments   &field_phi,
                  Field_Ell       &field_ell) const
  {

    using ExecPolicy = typename Kripke::Arch::Policy_LTimes<AL>::ExecPolicy;

    auto sdom_al = getSdomAL(al, sdom_id);
 
    // Get dimensioning
    int num_directions = set_dir.size(sdom_id);
    int num_groups =     set_group.size(sdom_id);
    int num_moments =    set_moment.size(sdom_id);
    int num_zones =      set_zone.size(sdom_id);

//#ifdef KRIPKE_USE_CHAI
//#ifdef KRIPKE_USE_CUDA
//    // Move data to GPU
//    sdom_al.moveHtoD(field_psi);
//    sdom_al.moveHtoD(field_phi);
//    sdom_al.moveHtoD(field_ell);
//    printf( "RCC try move H to D\n" );
//#endif
//#endif

#ifdef KRIPKE_USE_CHAI
#ifdef KRIPKE_USE_CUDA
    // Move data to GPU
    sdom_al.moveHtoD(field_psi);
    sdom_al.moveHtoD(field_phi);
    sdom_al.moveHtoD(field_ell);
    printf( "RCC AGAIN try move H to D\n" );

    // Get pointers
    auto psi = sdom_al.getDeviceView(field_psi);
    auto phi = sdom_al.getDeviceView(field_phi);
    auto ell = sdom_al.getDeviceView(field_ell);
    printf( "RCC device views\n" );
#endif

#else
    // Get pointers
    auto psi = sdom_al.getView(field_psi);
    auto phi = sdom_al.getView(field_phi);
    auto ell = sdom_al.getView(field_ell);

#endif

    // Compute:  phi =  ell * psi
    RAJA::kernel<ExecPolicy>(
        camp::make_tuple(
            RAJA::TypedRangeSegment<Moment>(0, num_moments),
            RAJA::TypedRangeSegment<Direction>(0, num_directions),
            RAJA::TypedRangeSegment<Group>(0, num_groups),
            RAJA::TypedRangeSegment<Zone>(0, num_zones) ),
        KRIPKE_LAMBDA (Moment nm, Direction d, Group g, Zone z) {

           phi(nm,g,z) += ell(nm, d) * psi(d, g, z);
           //ell(nm,d) *= 2;

        }
    );

#ifdef KRIPKE_USE_CHAI
#ifdef KRIPKE_USE_CUDA
    // Move data to CPU
    sdom_al.moveDtoH(field_psi);
    sdom_al.moveDtoH(field_phi);
    sdom_al.moveDtoH(field_ell);
#endif
#endif


  }

};

const std::string LTimesSdom::KernelName = "LTimes";






void Kripke::Kernel::LTimes(Kripke::Core::DataStore &data_store)
{
  KRIPKE_TIMER(data_store, LTimes);

  ArchLayoutV al_v = data_store.getVariable<ArchLayout>("al").al_v;

  Set const &set_dir    = data_store.getVariable<Set>("Set/Direction");
  Set const &set_group  = data_store.getVariable<Set>("Set/Group");
  Set const &set_zone   = data_store.getVariable<Set>("Set/Zone");
  Set const &set_moment = data_store.getVariable<Set>("Set/Moment");

  auto &field_psi =       data_store.getVariable<Field_Flux>("psi");
  auto &field_phi =       data_store.getVariable<Field_Moments>("phi");
  auto &field_ell =       data_store.getVariable<Field_Ell>("ell");

  // Loop over Subdomains
  for (Kripke::SdomId sdom_id : field_psi.getWorkList()){


    Kripke::dispatch(al_v, LTimesSdom{}, sdom_id,
                     set_dir, set_group, set_zone, set_moment,
                     field_psi, field_phi, field_ell);


  }

}


