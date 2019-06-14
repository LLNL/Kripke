//
// Copyright (c) 2014-19, Lawrence Livermore National Security, LLC
// and Kripke project contributors. See the COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
//

#ifndef KRIPKE_STEADYSTATESOLVER_H__
#define KRIPKE_STEADYSTATESOLVER_H__


#include <Kripke/Core/DataStore.h>

namespace Kripke {

  class DataStore;

  int SteadyStateSolver(Kripke::Core::DataStore &data_store, size_t max_iter, bool block_jacobi);


} // namespace

#endif

