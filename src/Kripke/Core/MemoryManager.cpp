//
// Copyright (c) 2014-25, Lawrence Livermore National Security, LLC
// and Kripke project contributors. See the Kripke/COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
//

#include <Kripke.h>
#include <Kripke/Core/MemoryManager.h>

#ifdef KRIPKE_USE_CHAI
#define DEBUG
#include <umpire/Umpire.hpp>
#include <umpire/strategy/QuickPool.hpp>
#undef DEBUG
#endif

using namespace Kripke;
using namespace Kripke::Core;

MemoryManager::MemoryManager(int device_pool_size) : device_pool_size(device_pool_size) {
#ifdef KRIPKE_USE_CHAI
  auto &rm = umpire::ResourceManager::getInstance();
  const char * allocator_name = "KRIPKE_DEVICE_POOL";
  size_t umpire_device_pool_size = ((size_t) device_pool_size) * 1024 * 1024 * 1024;
  size_t umpire_dev_block_size = 512;
  auto device_pool_allocator = rm.makeAllocator<umpire::strategy::QuickPool>(allocator_name, rm.getAllocator("DEVICE"), umpire_device_pool_size, umpire_dev_block_size);

  // Force the pool to materialize during initialization instead of on first use.
  void *tmp = device_pool_allocator.allocate(100*sizeof(int));
  device_pool_allocator.deallocate(tmp);
#endif // KRIPKE_USE_CHAI
}

double MemoryManager::getDeviceMemoryPoolSize() {
#ifdef KRIPKE_USE_CHAI
  return (double) device_pool_size;
#else
      return 0.0;
#endif
}

double MemoryManager::getDeviceMemoryHighWatermark() {
#ifdef KRIPKE_USE_CHAI
  auto device_allocator = getDeviceAllocator();
  return ((double) device_allocator.getHighWatermark()) / (1024 * 1024 * 1024);
#else
  return 0.0;
#endif
}

#ifdef KRIPKE_USE_CHAI
umpire::Allocator MemoryManager::getHostAllocator() {
  auto &rm = umpire::ResourceManager::getInstance();
  return rm.getAllocator("HOST");
}

umpire::Allocator MemoryManager::getDeviceAllocator() {
  auto &rm = umpire::ResourceManager::getInstance();
  return rm.getAllocator("KRIPKE_DEVICE_POOL");
}

void MemoryManager::copy(void *dst, void const *src, size_t bytes) {
  if(bytes == 0){
    return;
  }

  auto &rm = umpire::ResourceManager::getInstance();
  rm.copy(dst, src, bytes);
}
#endif
