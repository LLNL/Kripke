//
// Copyright (c) 2014-19, Lawrence Livermore National Security, LLC
// and Kripke project contributors. See the COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
//

#ifndef KRIPKE_SWEEPSOLVER_H__
#define KRIPKE_SWEEPSOLVER_H__

#include <Kripke.h>
#include <Kripke/Core/DataStore.h>
#include <Kripke/VarTypes.h>
#include <vector>

namespace Kripke {

  class DataStore;

  void SweepSolver (Kripke::Core::DataStore &data_store,
                    std::vector<SdomId> subdomain_list,
                    bool block_jacobi);



} // namespace

#endif

