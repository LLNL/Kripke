//
// Copyright (c) 2014-19, Lawrence Livermore National Security, LLC
// and Kripke project contributors. See the COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
//


#ifndef KRIPKE_H__
#define KRIPKE_H__

#include <KripkeConfig.h>

#include <RAJA/RAJA.hpp>

#include <string>
#include <vector>
#include <stdio.h>
#include <cmath>
#include <strings.h>
#include <stdlib.h>

// Make sure that there's openmp support, otherwise error out
#ifdef KRIPKE_USE_OPENMP
#ifndef _OPENMP
#error "OpenMP selected for build, but OpenMP is not available"
#endif
#endif

#ifdef KRIPKE_USE_MPI
#include <mpi.h>
#endif

// Forward Decl
struct Grid_Data;

#define KRESTRICT __restrict__

#ifdef KRIPKE_USE_MPI
#define KRIPKE_ABORT(...) \
  printf(__VA_ARGS__); \
  MPI_Abort(MPI_COMM_WORLD, 1);
#else
#define KRIPKE_ABORT(...) \
  printf(__VA_ARGS__); \
  exit(1);
#endif


#define KRIPKE_ASSERT(EXPR, ...) \
  if(!(EXPR)){\
    KRIPKE_ABORT("Assertion Failed: " __VA_ARGS__); \
  }
   

#define KRIPKE_LAMBDA [=] RAJA_HOST_DEVICE

namespace Kripke {

  /**
   * Index used to specify a local subdomain
   */
  RAJA_INDEX_VALUE(SdomId, "SdomId");


  /**
   * Index used to specify a global subdomain
   */
  RAJA_INDEX_VALUE(GlobalSdomId, "GlobalSdomId");


}



/**
  Tags for which parallel algorithm to use.
*/
enum ParallelMethod {
  PMETHOD_SWEEP,
  PMETHOD_BJ
};



/**
 * Import RAJA types into Kripke::Arch to make defining policies a lot
 * cleaner
 */
namespace Kripke {
namespace Arch {

  using RAJA::loop_exec;
  using RAJA::seq_exec;
  using RAJA::simd_exec;
  using RAJA::seq_reduce;
  using RAJA::auto_atomic;
  using RAJA::seq_atomic;
  using RAJA::ArgList;
  using RAJA::KernelPolicy;
  using RAJA::statement::Collapse;
  using RAJA::statement::If;
  using RAJA::statement::Param;
  using RAJA::statement::Not;
  using RAJA::statement::For;
  using RAJA::statement::Hyperplane;
  using RAJA::statement::Lambda;
  using RAJA::statement::Tile;
  using RAJA::tile_fixed;

#ifdef KRIPKE_USE_OPENMP
  using RAJA::omp_parallel_collapse_exec;
  using RAJA::omp_parallel_for_exec;
  using RAJA::omp_reduce;
#endif

#ifdef KRIPKE_USE_CUDA
  using RAJA::cuda_exec;
  using RAJA::cuda_block_x_loop;
  using RAJA::cuda_block_y_loop;
  using RAJA::cuda_block_z_loop;
  using RAJA::cuda_thread_x_direct;
  using RAJA::cuda_thread_y_direct;
  using RAJA::cuda_thread_z_direct;
  using RAJA::cuda_thread_x_loop;
  using RAJA::cuda_thread_y_loop;
  using RAJA::cuda_thread_z_loop;
  using RAJA::cuda_reduce;
  using RAJA::cuda_atomic;
  using RAJA::statement::CudaKernel;
  using RAJA::statement::CudaKernelAsync;
  using RAJA::statement::CudaSyncThreads;
#endif

} // namespace Arch
} // namespace Kripke



#endif

