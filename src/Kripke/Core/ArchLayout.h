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

namespace Kripke {
namespace Core {


struct Layout_DGZ {};

using Layout_Default = Layout_DGZ;

struct Arch_Sequential {};
struct Arch_OpenMP {};
struct Arch_CUDA {};


template<typename ARCH, typename LAYOUT>
struct ArchLayout {
  using Arch = ARCH;
  using Layout = LAYOUT;
};





template<typename Layout, typename IdxTuple>
struct FieldLayout;

/*
 * Default layouts are all defined as the default layout with
 * stride-1 access on right-most index
 */
template<typename ... IdxTypes>
struct FieldLayout<Layout_Default,
                   camp::tuple<IdxTypes...>>{

  static constexpr size_t num_dims = sizeof...(IdxTypes);

  static constexpr ptrdiff_t stride_one_dim = (ptrdiff_t)(num_dims) - 1;

  static std::array<ptrdiff_t, num_dims> getPermutation(){
    RAJA::as_array<RAJA::MakePerm<num_dims>>::get();
  }


};



} } // namespace

#endif

