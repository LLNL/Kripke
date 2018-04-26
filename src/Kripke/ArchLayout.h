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

#ifndef KRIPKE_ARCHLAYOUT_H__
#define KRIPKE_ARCHLAYOUT_H__

#include <Kripke.h>

namespace Kripke {

struct ArchT_Sequential {};
struct ArchT_OpenMP {};
struct ArchT_CUDA {};

enum ArchV {
  ArchV_Sequential,
  ArchV_OpenMP,
  ArchV_CUDA
};

struct LayoutT_DGZ {};
struct LayoutT_DZG {};
struct LayoutT_GDZ {};
struct LayoutT_GZD {};
struct LayoutT_ZDG {};
struct LayoutT_ZGD {};

enum LayoutV {
  LayoutV_DGZ,
  LayoutV_DZG,
  LayoutV_GDZ,
  LayoutV_GZD,
  LayoutV_ZDG,
  LayoutV_ZGD
};

template<typename ARCH, typename LAYOUT>
struct ArchLayoutT {
  using arch_t = ARCH;
  using layout_t = LAYOUT;
};

struct ArchLayoutV {
  ArchV arch_v;
  LayoutV layout_v;
};



template<typename Function, typename ... Args>
RAJA_INLINE
void dispatchLayout(LayoutV layout_v, Function const &fcn, Args &&... args)
{
  switch(layout_v){
    case LayoutV_DGZ: fcn(LayoutT_DGZ{}, std::forward<Args>(args)...); break;
/*    case LayoutV_DZG: fcn(LayoutT_DZG{}, std::forward<Args>(args)...); break;
    case LayoutV_GDZ: fcn(LayoutT_GDZ{}, std::forward<Args>(args)...); break;
    case LayoutV_GZD: fcn(LayoutT_GZD{}, std::forward<Args>(args)...); break;
    case LayoutV_ZDG: fcn(LayoutT_ZDG{}, std::forward<Args>(args)...); break;
    case LayoutV_ZGD: fcn(LayoutT_ZGD{}, std::forward<Args>(args)...); break;
*/  default: break;
  } 
}

template<typename Function, typename ... Args>
RAJA_INLINE
void dispatchArch(ArchV arch_v, Function const &fcn, Args &&... args)
{
  switch(arch_v){
    case ArchV_Sequential: fcn(ArchT_Sequential{}, std::forward<Args>(args)...); break;
    case ArchV_OpenMP: fcn(ArchT_OpenMP{}, std::forward<Args>(args)...); break;
    case ArchV_CUDA: fcn(ArchT_CUDA{}, std::forward<Args>(args)...); break;
  }
}


template<typename arch_t>
struct DispatchHelper{

  template<typename layout_t, typename Function, typename ... Args>
  void operator()(layout_t, Function const &fcn, Args &&... args) const {
    using al_t = ArchLayoutT<arch_t, layout_t>;
    fcn(al_t{}, std::forward<Args>(args)...);
  }
};

template<typename Function, typename ... Args>
RAJA_INLINE
void dispatch(ArchLayoutV al_v, Function const &fcn, Args &&... args)
{
  dispatchArch(al_v.arch_v, [&](auto arch_t){
    DispatchHelper<decltype(arch_t)> helper;

    dispatchLayout(al_v.layout_v, helper, fcn, std::forward<Args>(args)...);
  });
}

} // namespace

#endif

