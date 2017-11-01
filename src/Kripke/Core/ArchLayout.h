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

#ifndef KRIPKE_CORE_ARCH_LAYOUT_H__
#define KRIPKE_CORE_ARCH_LAYOUT_H__

#include <Kripke.h>
#include <array>


namespace Kripke {
namespace Core {


struct Layout_DGZ {};

using Layout_Default = Layout_DGZ;

struct Arch_Sequential {};
struct Arch_OpenMP {};
struct Arch_CUDA {};


template<typename ARCH, typename LAYOUT_FAMILY>
struct ArchLayout {
  using Arch = ARCH;
  using LayoutFamily = LAYOUT_FAMILY;
};


using AL_Default = ArchLayout<Layout_DGZ, Arch_Sequential>;

/*
 * Default layout is canonical ordering.
 *
 * This class is specialized for fields that needs data layouts to change.
 */

template<typename LAYOUT_FAMILY, typename ... IndexTypes>
struct LayoutInfo {
  // Default stride-one-index is the right-most index
  constexpr static ptrdiff_t num_dims = sizeof...(IndexTypes);
  constexpr static ptrdiff_t stride_one_dim = ((ptrdiff_t)num_dims)-1;

  using LayoutFamily = LAYOUT_FAMILY;
  using Layout = RAJA::TypedLayout<RAJA::Index_type, camp::tuple<IndexTypes...>, stride_one_dim>;

  static std::array<RAJA::Index_type, num_dims> getPermutation(){
    return RAJA::as_array<RAJA::MakePerm<num_dims>>::get();
  }
};


template<typename LayoutFamily, typename ... IndexTypes>
using LayoutType = typename LayoutInfo<LayoutFamily, IndexTypes...>::Layout;

template<typename LayoutFamily, typename ElementType, typename ... IndexTypes>
using ViewType = RAJA::View<ElementType, LayoutType<LayoutFamily, IndexTypes...>>;



} } // namespace

#endif

