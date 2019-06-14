//
// Copyright (c) 2014-19, Lawrence Livermore National Security, LLC
// and Kripke project contributors. See the COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
//

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


void Kripke::Generate::generateEnergy(Kripke::Core::DataStore &data_store,
    InputVariables const &input_vars)
{

  PartitionSpace &pspace = data_store.getVariable<PartitionSpace>("pspace");

  // Create sets for energy discretization
  size_t ngrp_per_sdom = input_vars.num_groups /
                         pspace.getGlobalNumSubdomains(SPACE_P);

  std::vector<size_t> local_grps(pspace.getNumSubdomains(SPACE_P),
                                 ngrp_per_sdom);

  RangeSet *grp_set = new RangeSet(pspace, SPACE_P, local_grps);
  data_store.addVariable("Set/Group", grp_set);

  GlobalRangeSet *global_grp_set = new GlobalRangeSet(pspace, *grp_set);
  data_store.addVariable("Set/GlobalGroup", global_grp_set);


}


