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

#include <Kripke.h>
#include <Kripke/Arch/LTimes.h>
#include <Kripke/Kernel.h>
#include <Kripke/Timing.h>
#include <Kripke/VarTypes.h>

#include<utility>

using namespace Kripke;
using namespace Kripke::Core;


struct LTimesSdom {

  template<typename AL>
  RAJA_INLINE
  void operator()(AL, Kripke::Core::DataStore &,
                  Kripke::SdomId sdom_id,
                  Set const       &set_dir,
                  Set const       &set_group,
                  Set const       &set_zone,
                  Set const       &set_moment,
                  Field_Flux      &field_psi,
                  Field_Moments   &field_phi,
                  Field_Ell       &field_ell) const
  {
    // Get dimensioning
    int num_directions = set_dir.size(sdom_id);
    int num_groups =     set_group.size(sdom_id);
    int num_moments =    set_moment.size(sdom_id);
    int num_zones =      set_zone.size(sdom_id);

    // Get pointers
    auto psi = field_psi.getViewAL<AL>(sdom_id);
    auto phi = field_phi.getViewAL<AL>(sdom_id);
    auto ell = field_ell.getViewAL<AL>(sdom_id);

    // Compute:  phi =  ell * psi
		RAJA::kernel<Kripke::Arch::Policy_LTimes>(
        camp::make_tuple(
            RAJA::TypedRangeSegment<Moment>(0, num_moments),
            RAJA::TypedRangeSegment<Direction>(0, num_directions),
            RAJA::TypedRangeSegment<Group>(0, num_groups),
            RAJA::TypedRangeSegment<Zone>(0, num_zones) ),
        KRIPKE_LAMBDA (Moment nm, Direction d, Group g, Zone z) {

           phi(nm,g,z) += ell(nm, d) * psi(d, g, z);

        }
    );


  }

};








void Kripke::Kernel::LTimes(Kripke::Core::DataStore &data_store)
{
  KRIPKE_TIMER(data_store, LTimes);

  Set const &set_dir    = data_store.getVariable<Set>("Set/Direction");
  Set const &set_group  = data_store.getVariable<Set>("Set/Group");
  Set const &set_zone   = data_store.getVariable<Set>("Set/Zone");
  Set const &set_moment = data_store.getVariable<Set>("Set/Moment");

  auto &field_psi =       data_store.getVariable<Field_Flux>("psi");
  auto &field_phi =       data_store.getVariable<Field_Moments>("phi");
  auto &field_ell =       data_store.getVariable<Field_Ell>("ell");

  // Loop over Subdomains
  for (Kripke::SdomId sdom_id : field_psi.getWorkList()){


    Kripke::Kernel::dispatch(data_store, sdom_id, LTimesSdom{},
                             set_dir, set_group, set_zone, set_moment,
                             field_psi, field_phi, field_ell);


  }


}


