//
// Copyright (c) 2014-25, Lawrence Livermore National Security, LLC
// and Kripke project contributors. See the Kripke/COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
//

#ifndef KRIPKE_CORE_MEMORYMANAGER_H__
#define KRIPKE_CORE_MEMORYMANAGER_H__

#include <Kripke.h>

#ifdef KRIPKE_USE_CHAI
#include <umpire/Umpire.hpp>
#endif

namespace Kripke {
namespace Core {

class MemoryManager {
  protected:
    int device_pool_size;

  public:
    MemoryManager(int device_pool_size);
    double getDeviceMemoryPoolSize();
    double getDeviceMemoryHighWatermark();

#ifdef KRIPKE_USE_CHAI
    static umpire::Allocator getHostAllocator();
    static umpire::Allocator getDeviceAllocator();
    static void copy(void *dst, void const *src, size_t bytes);
#endif
};

} } // namespace

#endif
