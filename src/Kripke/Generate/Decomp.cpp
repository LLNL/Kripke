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

#include <Kripke/Generate.h>

#include <Kripke/Core/Comm.h>
#include <Kripke/Core/Field.h>
#include <Kripke/Kernel.h>
#include <Kripke/Core/PartitionSpace.h>
#include <Kripke/Core/Set.h>
#include <Kripke/Timing.h>
#include <Kripke/VarTypes.h>

using namespace Kripke;
using namespace Kripke::Core;


void Kripke::Generate::generateDecomp(Kripke::Core::DataStore &data_store,
    InputVariables const &input_vars)
{
  // Create a "Comm World"
  auto &comm = data_store.newVariable<Kripke::Core::Comm>("comm");

  // Create our partitioning over MPI
  auto &pspace = data_store.newVariable<PartitionSpace>("pspace",
      comm,
      1,
      1,
      input_vars.npx,
      input_vars.npy,
      input_vars.npz);

  // Create our local partition over subdomains
  pspace.setup_createSubdomains(
      input_vars.num_groupsets,
      input_vars.num_dirsets,
      input_vars.num_zonesets_dim[0],
      input_vars.num_zonesets_dim[1],
      input_vars.num_zonesets_dim[2]);

  // Create utility Sets and Fields that describe our global subdomain layout
  pspace.createSubdomainData(data_store);
  pspace.print();


}


