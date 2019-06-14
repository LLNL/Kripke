//
// Copyright (c) 2014-19, Lawrence Livermore National Security, LLC
// and Kripke project contributors. See the COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
//

#ifndef KRIPKE_ARCHLAYOUT_H__
#define KRIPKE_ARCHLAYOUT_H__

#include <Kripke.h>
#include <Kripke/Core/BaseVar.h>

#include <strings.h>

namespace Kripke {

struct ArchT_Sequential {};

#ifdef KRIPKE_USE_OPENMP
struct ArchT_OpenMP {};
#endif

#ifdef KRIPKE_USE_CUDA
struct ArchT_CUDA {};
#endif



enum ArchV {
  ArchV_Unknown = -1,
  ArchV_Sequential,

#ifdef KRIPKE_USE_OPENMP
  ArchV_OpenMP,
#endif

#ifdef KRIPKE_USE_CUDA
  ArchV_CUDA,
#endif

  ArchV_num_values
};


RAJA_INLINE
std::string archToString(ArchV av){
  switch(av){
    case ArchV_Sequential:    return "Sequential";

#ifdef KRIPKE_USE_OPENMP
    case ArchV_OpenMP:        return "OpenMP";
#endif

#ifdef KRIPKE_USE_CUDA
    case ArchV_CUDA:          return "CUDA";
#endif

    case ArchV_Unknown:
    case ArchV_num_values:
    default:                  return "unknown";
  }
}

RAJA_INLINE
ArchV stringToArch(std::string const &str){
  for(int av = 0;av < (int)ArchV_num_values;++ av){
    if(!strcasecmp(archToString((ArchV)av).c_str(), str.c_str())){
      return (ArchV)av;
    }
  }
  return ArchV_Unknown;
}

struct LayoutT_DGZ {};
struct LayoutT_DZG {};
struct LayoutT_GDZ {};
struct LayoutT_GZD {};
struct LayoutT_ZDG {};
struct LayoutT_ZGD {};

enum LayoutV {
  LayoutV_Unknown = -1,
  LayoutV_DGZ,
  LayoutV_DZG,
  LayoutV_GDZ,
  LayoutV_GZD,
  LayoutV_ZDG,
  LayoutV_ZGD,
  LayoutV_num_values
};

RAJA_INLINE
std::string layoutToString(LayoutV lv){
  switch(lv){
    case LayoutV_DGZ:    return "DGZ";
    case LayoutV_DZG:    return "DZG";
    case LayoutV_GDZ:    return "GDZ";
    case LayoutV_GZD:    return "GZD";
    case LayoutV_ZDG:    return "ZDG";
    case LayoutV_ZGD:    return "ZGD";
    case LayoutV_Unknown:
    case LayoutV_num_values:
    default:             return "unknown";
  }
}

RAJA_INLINE
LayoutV stringToLayout(std::string const &str){
  for(int lv = 0;lv < (int)LayoutV_num_values;++ lv){
    if(!strcasecmp(layoutToString((LayoutV)lv).c_str(), str.c_str())){
      return (LayoutV)lv;
    }
  }
  return LayoutV_Unknown;
}


template<typename ARCH, typename LAYOUT>
struct ArchLayoutT {
  using arch_t = ARCH;
  using layout_t = LAYOUT;
};

struct ArchLayoutV {
  ArchV arch_v;
  LayoutV layout_v;
};


class ArchLayout : public Kripke::Core::BaseVar {
public:
  ArchLayout() = default;
  virtual ~ArchLayout() = default;

  ArchLayoutV al_v;
};


template<typename Function, typename ... Args>
RAJA_INLINE
void dispatchLayout(LayoutV layout_v, Function const &fcn, Args &&... args)
{
  switch(layout_v){
    case LayoutV_DGZ: fcn(LayoutT_DGZ{}, std::forward<Args>(args)...); break;
    case LayoutV_DZG: fcn(LayoutT_DZG{}, std::forward<Args>(args)...); break;
    case LayoutV_GDZ: fcn(LayoutT_GDZ{}, std::forward<Args>(args)...); break;
    case LayoutV_GZD: fcn(LayoutT_GZD{}, std::forward<Args>(args)...); break;
    case LayoutV_ZDG: fcn(LayoutT_ZDG{}, std::forward<Args>(args)...); break;
    case LayoutV_ZGD: fcn(LayoutT_ZGD{}, std::forward<Args>(args)...); break;
    default: KRIPKE_ABORT("Unknown layout_v=%d\n", (int)layout_v); break;
  } 
}

template<typename Function, typename ... Args>
RAJA_INLINE
void dispatchArch(ArchV arch_v, Function const &fcn, Args &&... args)
{
  switch(arch_v){
    case ArchV_Sequential: fcn(ArchT_Sequential{}, std::forward<Args>(args)...); break;
#ifdef KRIPKE_USE_OPENMP
    case ArchV_OpenMP: fcn(ArchT_OpenMP{}, std::forward<Args>(args)...); break;
#endif 

#ifdef KRIPKE_USE_CUDA
    case ArchV_CUDA: fcn(ArchT_CUDA{}, std::forward<Args>(args)...); break;
#endif
    default: KRIPKE_ABORT("Unknown arch_v=%d\n", (int)arch_v); break;
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

