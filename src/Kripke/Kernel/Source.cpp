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

#include <Kripke.h>
#include <Kripke/Arch/Source.h>
#include <Kripke/Timing.h>
#include <Kripke/VarTypes.h>

/**
 * Add an isotropic source, with flux of 1, to every zone with Region 1
 * (or material 0).
 *
 * Since it's isotropic, we're just adding this to nm=0.
 */
void Kripke::Kernel::source(Kripke::Core::DataStore &data_store)
{
  KRIPKE_TIMER(data_store, Source);

  auto &set_group   = data_store.getVariable<Kripke::Core::Set>("Set/Group");
  auto &set_mixelem = data_store.getVariable<Kripke::Core::Set>("Set/MixElem");

  auto &field_phi_out = data_store.getVariable<Kripke::Field_Moments>("phi_out");

  auto &field_mixed_to_zone     = data_store.getVariable<Field_MixElem2Zone>("mixelem_to_zone");
  auto &field_mixed_to_material = data_store.getVariable<Field_MixElem2Material>("mixelem_to_material");
  auto &field_mixed_to_fraction = data_store.getVariable<Field_MixElem2Double>("mixelem_to_fraction");

  // Source term is isotropic
  Moment nm{0};

  // Loop over zoneset subdomains
  for(auto sdom_id : field_phi_out.getWorkList()){

    auto mixelem_to_zone     = field_mixed_to_zone.getView(sdom_id);
    auto mixelem_to_material = field_mixed_to_material.getView(sdom_id);
    auto mixelem_to_fraction = field_mixed_to_fraction.getView(sdom_id);

    auto phi_out = field_phi_out.getView(sdom_id);

    int num_mixed  = set_mixelem.size(sdom_id);
    int num_groups = set_group.size(sdom_id);


    // Compute:  phi =  ell * psi
    RAJA::nested::forall(Kripke::Arch::Policy_Source{},
        RAJA::util::make_tuple(
            RAJA::RangeSegment(0, num_groups),
            RAJA::RangeSegment(0, num_mixed) ),
        KRIPKE_LAMBDA (Group g, MixElem mix) {

            Material material = mixelem_to_material(mix);

            if(material == 2){
              Zone z = mixelem_to_zone(mix);
              double fraction = mixelem_to_fraction(mix);

              phi_out(nm, g, z) += 1.0 * fraction;
            }

        }
    );

  }
}

