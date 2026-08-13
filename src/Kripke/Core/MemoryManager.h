//
// Copyright (c) 2014-25, Lawrence Livermore National Security, LLC
// and Kripke project contributors. See the Kripke/COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
//

#ifndef KRIPKE_CORE_MEMORYMANAGER_H__
#define KRIPKE_CORE_MEMORYMANAGER_H__

#include <Kripke.h>
#include <cstddef>

namespace Kripke {
namespace Core {

class MemoryManager {
  protected:
    size_t requested_device_pool_size;

  public:
    MemoryManager(size_t requested_device_pool_size);
    double getDeviceMemoryPoolSize();
    double getDeviceMemoryHighWatermark();
};

} } // namespace

#endif
