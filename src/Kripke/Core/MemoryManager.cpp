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
#include <chai/ManagedArray.hpp>
#undef DEBUG
#endif

using namespace Kripke;
using namespace Kripke::Core;

MemoryManager::MemoryManager(int device_pool_size) : device_pool_size(device_pool_size) {
#ifdef KRIPKE_USE_CHAI && (defined(KRIPKE_USE_CUDA) || defined(KRIPKE_USE_HIP))
  auto &rm = umpire::ResourceManager::getInstance();
  const char * allocator_name = "KRIPKE_DEVICE_POOL";
  size_t umpire_device_pool_size = ((size_t) device_pool_size) * 1024 * 1024 * 1024;
  size_t umpire_dev_block_size = 512;
  auto device_pool_allocator = rm.makeAllocator<umpire::strategy::QuickPool>(allocator_name, rm.getAllocator("DEVICE"), umpire_device_pool_size, umpire_dev_block_size);
  auto chai_resource_manager = chai::ArrayManager::getInstance();
  chai_resource_manager->setAllocator(chai::GPU, device_pool_allocator);
  // force allocation of GPU memory pool
  auto tmp = new chai::ManagedArray<int>(100, chai::GPU);
  tmp->free(chai::GPU);
  delete tmp;
#endif // KRIPKE_USE_CHAI
}

double MemoryManager::getDeviceMemoryPoolSize() {
#ifdef KRIPKE_USE_CHAI && (defined(KRIPKE_USE_CUDA) || defined(KRIPKE_USE_HIP))
  return (double) device_pool_size;
#else
      return 0.0;
#endif
}

double MemoryManager::getDeviceMemoryHighWatermark() {
#ifdef KRIPKE_USE_CHAI && (defined(KRIPKE_USE_CUDA) || defined(KRIPKE_USE_HIP))
  auto chai_resource_manager = chai::ArrayManager::getInstance();
  auto device_allocator = chai_resource_manager->getAllocator(chai::GPU);
  return ((double) device_allocator.getHighWatermark()) / (1024 * 1024 * 1024);
#else
  return 0.0;
#endif
}
