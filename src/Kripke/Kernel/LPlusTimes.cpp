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
#include <Kripke/Arch/LPlusTimes.h>
#include <Kripke/Kernel.h>
#include <Kripke/Timing.h>
#include <Kripke/VarTypes.h>

using namespace Kripke;
using namespace Kripke::Core;

struct LPlusTimesSdom {

  template<typename AL>
  void operator()(AL, Kripke::Core::DataStore &,
                  Kripke::SdomId sdom_id,
                  Set const       &set_dir,
                  Set const       &set_group,
                  Set const       &set_zone,
                  Set const       &set_moment,
                  Field_Moments   &field_phi_out,
                  Field_Flux      &field_rhs,
                  Field_EllPlus   &field_ell_plus) const
  {

    // Get dimensioning
    int num_directions = set_dir.size(sdom_id);
    int num_groups =     set_group.size(sdom_id);
    int num_moments =    set_moment.size(sdom_id);
    int num_zones =      set_zone.size(sdom_id);

    // Get views
    auto phi_out =  field_phi_out.getViewAL<AL>(sdom_id);
    auto rhs =      field_rhs.getViewAL<AL>(sdom_id);
    auto ell_plus = field_ell_plus.getViewAL<AL>(sdom_id);

    // Compute:  rhs =  ell_plus * phi_out
    RAJA::kernel<Kripke::Arch::Policy_LPlusTimes>(
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

  // Loop over Subdomains
  for (Kripke::SdomId sdom_id : field_rhs.getWorkList()){

    Kripke::Kernel::dispatch(data_store, sdom_id, LPlusTimesSdom{},
                             set_dir, set_group, set_zone, set_moment,
                             field_phi_out, field_rhs, field_ell_plus);

  }


}
