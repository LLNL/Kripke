/*
 * NOTICE
 *
 * This work was produced at the Lawrence Livermore National Laboratory (LLNL)
 * under contract no. DE-AC-52-07NA27344 (Contract 44) between the U.S.
 * Department of Energy (DOE) and Lawrence Livermore National Security, LLC
 * (LLNS) for the operation of LLNL. The rights of the Federal Government are
 * reserved under Contract 44.
 *
 * DISCLAIMER
 *
 * This work was prepared as an account of work sponsored by an agency of the
 * United States Government. Neither the United States Government nor Lawrence
 * Livermore National Security, LLC nor any of their employees, makes any
 * warranty, express or implied, or assumes any liability or responsibility
 * for the accuracy, completeness, or usefulness of any information, apparatus,
 * product, or process disclosed, or represents that its use would not infringe
 * privately-owned rights. Reference herein to any specific commercial products,
 * process, or service by trade name, trademark, manufacturer or otherwise does
 * not necessarily constitute or imply its endorsement, recommendation, or
 * favoring by the United States Government or Lawrence Livermore National
 * Security, LLC. The views and opinions of authors expressed herein do not
 * necessarily state or reflect those of the United States Government or
 * Lawrence Livermore National Security, LLC, and shall not be used for
 * advertising or product endorsement purposes.
 *
 * NOTIFICATION OF COMMERCIAL USE
 *
 * Commercialization of this product is prohibited without notifying the
 * Department of Energy (DOE) or Lawrence Livermore National Security.
 */

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
  using RAJA::atomic::auto_atomic;
  using RAJA::atomic::seq_atomic;
  using RAJA::ArgList;
  using RAJA::KernelPolicy;
  using RAJA::statement::Collapse;
  using RAJA::statement::If;
  using RAJA::statement::Param;
  using RAJA::statement::Not;
  using RAJA::statement::For;
  using RAJA::statement::Hyperplane;
  using RAJA::statement::Lambda;
  using RAJA::statement::SetShmemWindow;
  using RAJA::statement::Tile;
  using RAJA::statement::tile_fixed;

#ifdef KRIPKE_USE_OPENMP
  using RAJA::omp_parallel_collapse_exec;
  using RAJA::omp_parallel_for_exec;
  using RAJA::omp_reduce;
#endif

} // namespace Arch
} // namespace Kripke



#endif

