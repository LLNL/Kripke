//
// Copyright (c) 2014-25, Lawrence Livermore National Security, LLC
// and Kripke project contributors. See the Kripke/COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
//

#include <Kripke.h>
#include <Kripke/Core/MemoryManager.h>

using namespace Kripke;
using namespace Kripke::Core;

MemoryManager::MemoryManager(int device_pool_size) : device_pool_size(device_pool_size) {
#if defined(KRIPKE_USE_UMPIRE) && (defined(KRIPKE_USE_CUDA) || defined(KRIPKE_USE_HIP))
  auto &rm = umpire::ResourceManager::getInstance();
  const char * allocator_name = "KRIPKE_DEVICE_POOL";
  size_t umpire_device_pool_size = ((size_t) device_pool_size) * 1024 * 1024 * 1024;
  size_t umpire_dev_block_size = 512;
  auto device_pool_allocator = rm.makeAllocator<umpire::strategy::QuickPool>(allocator_name, rm.getAllocator("DEVICE"), umpire_device_pool_size, umpire_dev_block_size);
  // Force allocation of GPU memory pool
  void *tmp = device_pool_allocator.allocate(100*sizeof(int));
  device_pool_allocator.deallocate(tmp);
#if defined(KRIPKE_USE_CHAI)
  // Set CHAI device memory pool allocator
  auto chai_resource_manager = chai::ArrayManager::getInstance();
  chai_resource_manager->setAllocator(Kripke::GPU, device_pool_allocator);
#endif // KRIPKE_USE_CHAI
#endif // KRIPKE_USE_UMPIRE
}

double MemoryManager::getDeviceMemoryPoolSize() {
#if defined(KRIPKE_USE_UMPIRE) && (defined(KRIPKE_USE_CUDA) || defined(KRIPKE_USE_HIP))
  return (double) device_pool_size;
#else
      return 0.0;
#endif
}

double MemoryManager::getDeviceMemoryHighWatermark() {
#if defined(KRIPKE_USE_UMPIRE) && (defined(KRIPKE_USE_CUDA) || defined(KRIPKE_USE_HIP))
  auto device_allocator = getDeviceAllocator();
  return ((double) device_allocator.getHighWatermark()) / (1024 * 1024 * 1024);
#else
  return 0.0;
#endif
}

#if defined(KRIPKE_USE_UMPIRE)
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
  rm.copy(dst, const_cast<void *>(src), bytes);
}
#endif
