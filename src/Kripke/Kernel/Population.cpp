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

#include <Kripke/Comm.h>
#include <Kripke/Grid.h>
#include <Kripke/SubTVec.h>
#include <Kripke/Timing.h>
#include <Kripke/VarTypes.h>

/**
 * Returns the integral of Psi over all phase-space, to look at convergence
 */
double Kripke::Kernel::population(Kripke::DataStore &data_store)
{
  KRIPKE_TIMER(data_store, Population);

  Set const &set_dir    = data_store.getVariable<Set>("Set/Direction");
  Set const &set_group  = data_store.getVariable<Kripke::Set>("Set/Group");
  Set const &set_zone   = data_store.getVariable<Set>("Set/Zone");

  auto &field_psi = data_store.getVariable<Kripke::Field_Flux>("psi");
  auto &field_w = data_store.getVariable<Field_Direction2Double>("quadrature/w");
  auto &field_volume = data_store.getVariable<Field_Zone2Double>("volume");

  // sum up particles for psi and rhs
  double part = 0.0;
  for (Kripke::SdomId sdom_id : field_psi.getWorkList()){

    int num_directions = set_dir.size(sdom_id);
    int num_groups =     set_group.size(sdom_id);
    int num_zones =      set_zone.size(sdom_id);

    auto psi = field_psi.getView(sdom_id);
    auto w = field_w.getView(sdom_id);
    auto volume = field_volume.getView(sdom_id);

    for(Zone z{0};z < num_zones;++ z){
      for(Direction d{0};d < num_directions;++ d){
        for(Group g{0};g < num_groups;++ g){

          part += w(d) * psi(d,g,z) * volume(z);

        }
      }
    }
  }

  // reduce
  auto const &comm = data_store.getVariable<Kripke::Comm>("comm");
  return comm.allReduceSumDouble(part);
}

