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

#ifndef KRIPKE_ARCH_POPULATION
#define KRIPKE_ARCH_POPULATION

#include <Kripke.h>
#include <Kripke/VarTypes.h>

namespace Kripke {
namespace Arch {

#ifdef KRIPKE_ARCH_SEQUENTIAL

using Reduce_Population = RAJA::seq_reduce;

using Policy_Population =
    RAJA::KernelPolicy<
      RAJA::statement::For<0, RAJA::loop_exec,
        RAJA::statement::For<1, RAJA::loop_exec,
          RAJA::statement::For<2, RAJA::loop_exec,
            RAJA::statement::Lambda<0>
          >
        >
      >
    >;
#endif

#ifdef KRIPKE_ARCH_OPENMP

using Reduce_Population = RAJA::omp_reduce;

using Policy_Population =
    RAJA::KernelPolicy<
      RAJA::statement::For<1, RAJA::omp_parallel_for_exec,
        RAJA::statement::For<0, RAJA::loop_exec,
          RAJA::statement::For<2, RAJA::loop_exec,
            RAJA::statement::Lambda<0>
          >
        >
      >
    >;
#endif

}
}

#endif
