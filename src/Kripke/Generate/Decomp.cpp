//
// Copyright (c) 2014-19, Lawrence Livermore National Security, LLC
// and Kripke project contributors. See the COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
//

#include <Kripke/Generate.h>

#include <Kripke/Core/Comm.h>
#include <Kripke/Core/Field.h>
#include <Kripke/ArchLayout.h>
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

  // Create our ArchLayout object to describe how we are going to 
  // execute, and what data layouts we want
  auto &al_var = data_store.newVariable<ArchLayout>("al");
  al_var.al_v = input_vars.al_v;

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


