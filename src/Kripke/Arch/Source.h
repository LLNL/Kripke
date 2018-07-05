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

#ifndef KRIPKE_ARCH_SOURCE
#define KRIPKE_ARCH_SOURCE

#include <Kripke.h>
#include <Kripke/VarTypes.h>

namespace Kripke {
namespace Arch {

template<typename AL>
struct Policy_Source;

template<typename A>
struct Policy_Source<ArchLayoutT<A, LayoutT_GDZ>> :
  Policy_Source<ArchLayoutT<A, LayoutT_DGZ>>{};

template<typename A>
struct Policy_Source<ArchLayoutT<A, LayoutT_GZD>> :
  Policy_Source<ArchLayoutT<A, LayoutT_DGZ>>{};

template<typename A>
struct Policy_Source<ArchLayoutT<A, LayoutT_ZDG>> :
  Policy_Source<ArchLayoutT<A, LayoutT_DZG>>{};

template<typename A>
struct Policy_Source<ArchLayoutT<A, LayoutT_ZGD>> :
  Policy_Source<ArchLayoutT<A, LayoutT_DZG>>{};

template<>
struct Policy_Source<ArchLayoutT<ArchT_Sequential, LayoutT_DGZ>> {
  using ExecPolicy = 
    KernelPolicy<
      Collapse<loop_exec, ArgList<0,1>, // Group, MixElem
        Lambda<0>
      >
    >;
};

template<>
struct Policy_Source<ArchLayoutT<ArchT_Sequential, LayoutT_DZG>> {
  using ExecPolicy = 
    KernelPolicy<
      Collapse<loop_exec, ArgList<0,1>, // MixElem, Group
        Lambda<0>
      >
    >;
};




#ifdef KRIPKE_USE_OPENMP
template<>
struct Policy_Source<ArchLayoutT<ArchT_OpenMP, LayoutT_DGZ>> {
  using ExecPolicy =
    KernelPolicy<
      Collapse<omp_parallel_collapse_exec, ArgList<0,1>, // Group, MixElem
        Lambda<0>
      >
    >;
};

template<>
struct Policy_Source<ArchLayoutT<ArchT_OpenMP, LayoutT_DZG>> {
  using ExecPolicy =
    KernelPolicy<
      Collapse<omp_parallel_collapse_exec, ArgList<1,0>, // MixElem, Group
        Lambda<0>
      >
    >;
};
#endif // KRIPKE_USE_OPENMP


#ifdef KRIPKE_USE_CUDA
template<>
struct Policy_Source<ArchLayoutT<ArchT_CUDA, LayoutT_DGZ>> {
  using ExecPolicy =
    KernelPolicy<
      Collapse<loop_exec, ArgList<0,1>, // Group, MixElem
        Lambda<0>
      >
    >;
};

template<>
struct Policy_Source<ArchLayoutT<ArchT_CUDA, LayoutT_DZG>> {
  using ExecPolicy =
    KernelPolicy<
      Collapse<loop_exec, ArgList<0,1>, // MixElem, Group
        Lambda<0>
      >
    >;
};
#endif // KRIPKE_USE_CUDA


}
}

#endif
