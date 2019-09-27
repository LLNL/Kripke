//
// Copyright (c) 2014-19, Lawrence Livermore National Security, LLC
// and Kripke project contributors. See the COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
//

#ifndef KRIPKE_ARCH_LTIMES
#define KRIPKE_ARCH_LTIMES

#include <Kripke.h>
#include <Kripke/VarTypes.h>

namespace Kripke {
namespace Arch {

template<typename AL>
struct Policy_LTimes;

template<>
struct Policy_LTimes<ArchLayoutT<ArchT_Sequential, LayoutT_DGZ>> {
  using ExecPolicy = 
    KernelPolicy<
      For<0, loop_exec, // moment
        For<1, loop_exec, // direction
          For<2, loop_exec, // group
            For<3, loop_exec, // zone
              Lambda<0>
            >
          >
        >
      >
    >;
};

template<>
struct Policy_LTimes<ArchLayoutT<ArchT_Sequential, LayoutT_DZG>> {
  using ExecPolicy = 
    KernelPolicy<
      For<0, loop_exec, // moment
        For<1, loop_exec, // direction
          For<3, loop_exec, // zone
            For<2, loop_exec, // group
              Lambda<0>
            >
          >
        >
      >
    >;
};

template<>
struct Policy_LTimes<ArchLayoutT<ArchT_Sequential, LayoutT_GDZ>> {
  using ExecPolicy = 
    KernelPolicy<
      For<2, loop_exec, // group
        For<0, loop_exec, // moment
          For<1, loop_exec, // direction
            For<3, loop_exec, // zone
              Lambda<0>
            >
          >
        >
      >
    >;
};

template<>
struct Policy_LTimes<ArchLayoutT<ArchT_Sequential, LayoutT_GZD>> {
  using ExecPolicy = 
    KernelPolicy<
      For<2, loop_exec, // group
        For<3, loop_exec, // zone
          For<0, loop_exec, // moment
            For<1, loop_exec, // direction
              Lambda<0>
            >
          >
        >
      >
    >;
};

template<>
struct Policy_LTimes<ArchLayoutT<ArchT_Sequential, LayoutT_ZDG>> {
  using ExecPolicy = 
    KernelPolicy<
      For<3, loop_exec, // zone
        For<0, loop_exec, // moment
          For<1, loop_exec, // direction
            For<2, loop_exec, // group
              Lambda<0>
            >
          >
        >
      >
    >;
};

template<>
struct Policy_LTimes<ArchLayoutT<ArchT_Sequential, LayoutT_ZGD>> {
  using ExecPolicy =
    KernelPolicy<
      For<3, loop_exec, // zone
        For<2, loop_exec, // group
          For<0, loop_exec, // moment
            For<1, loop_exec, // direction
              Lambda<0>
            >
          >
        >
      >
    >;
};




#ifdef KRIPKE_USE_OPENMP
template<>
struct Policy_LTimes<ArchLayoutT<ArchT_OpenMP, LayoutT_DGZ>> {
  using ExecPolicy =
    KernelPolicy<
      Collapse<omp_parallel_collapse_exec, ArgList<0,2>, // Moment Group
        For<1, loop_exec, // Direction
          For<3, loop_exec, // Zone
            Lambda<0>
          >
        >
      >
    >;
};

template<>
struct Policy_LTimes<ArchLayoutT<ArchT_OpenMP, LayoutT_DZG>> {
  using ExecPolicy =
    KernelPolicy<
      Collapse<omp_parallel_collapse_exec, ArgList<0,3>, // Moment Zone
        For<1, loop_exec, // Direction
          For<2, loop_exec, // Group
            Lambda<0>
          >
        >
      >
    >;
};

template<>
struct Policy_LTimes<ArchLayoutT<ArchT_OpenMP, LayoutT_GDZ>> {
  using ExecPolicy =
    KernelPolicy<
      Collapse<omp_parallel_collapse_exec, ArgList<2,0>, // Group Moment
        For<1, loop_exec, // Direction
          For<3, loop_exec, // Zone
            Lambda<0>
          >
        >
      >
    >;
};

template<>
struct Policy_LTimes<ArchLayoutT<ArchT_OpenMP, LayoutT_GZD>> {
  using ExecPolicy =
    KernelPolicy<
      Collapse<omp_parallel_collapse_exec, ArgList<2,3,0>, // Group Zone Moment
        For<1, loop_exec, // Direection
          Lambda<0>
        >
      >
    >;
};

template<>
struct Policy_LTimes<ArchLayoutT<ArchT_OpenMP, LayoutT_ZDG>> {
  using ExecPolicy =
    KernelPolicy<
      Collapse<omp_parallel_collapse_exec, ArgList<3,0>, // Zone Moment
        For<1, loop_exec, // Direction
          For<2, loop_exec, // Group
            Lambda<0>
          >
        >
      >
    >;
};

template<>
struct Policy_LTimes<ArchLayoutT<ArchT_OpenMP, LayoutT_ZGD>> {
  using ExecPolicy =
    KernelPolicy<
      Collapse<omp_parallel_collapse_exec, ArgList<3,2>, // Zone Group
        For<0, loop_exec, // Moment
          For<1, loop_exec, // Direction
            Lambda<0>
          >
        >
      >
    >;
};
#endif // KRIPKE_USE_OPENMP



#ifdef KRIPKE_USE_CUDA
template<>
struct Policy_LTimes<ArchLayoutT<ArchT_CUDA, LayoutT_DGZ>> {
  using ExecPolicy =
    KernelPolicy<
      CudaKernel<
        For<2, cuda_block_x_loop, // group
          For<0, cuda_block_y_loop, // moment
            For<3, cuda_thread_x_loop, // zone
              For<1, seq_exec, // direction
                Lambda<0>
              >
            >
          >
        >
      >
    >;
};

template<>
struct Policy_LTimes<ArchLayoutT<ArchT_CUDA, LayoutT_DZG>> {
    using ExecPolicy =
      KernelPolicy<
        CudaKernel<
          For<2, cuda_block_x_loop, // group
            For<0, cuda_block_y_loop, // moment
              For<3, cuda_thread_x_loop, // zone
                For<1, seq_exec, // direction
                  Lambda<0>
                >
              >
            >
          >
        >
      >;
};

template<>
struct Policy_LTimes<ArchLayoutT<ArchT_CUDA, LayoutT_GDZ>> {
    using ExecPolicy =
      KernelPolicy<
        CudaKernel<
          For<2, cuda_block_x_loop, // group
            For<0, cuda_block_y_loop, // moment
              For<3, cuda_thread_x_loop, // zone
                For<1, seq_exec, // direction
                  Lambda<0>
                >
              >
            >
          >
        >
      >;
};

template<>
struct Policy_LTimes<ArchLayoutT<ArchT_CUDA, LayoutT_GZD>> {
    using ExecPolicy =
      KernelPolicy<
        CudaKernel<
          For<2, cuda_block_x_loop, // group
            For<0, cuda_block_y_loop, // moment
              For<3, cuda_thread_x_loop, // zone
                For<1, seq_exec, // direction
                  Lambda<0>
                >
              >
            >
          >
        >
      >;
};

template<>
struct Policy_LTimes<ArchLayoutT<ArchT_CUDA, LayoutT_ZDG>> {
    using ExecPolicy =
      KernelPolicy<
        CudaKernel<
          For<2, cuda_block_x_loop, // group
            For<0, cuda_block_y_loop, // moment
              For<3, cuda_thread_x_loop, // zone
                For<1, seq_exec, // direction
                  Lambda<0>
                >
              >
            >
          >
        >
      >;
};

template<>
struct Policy_LTimes<ArchLayoutT<ArchT_CUDA, LayoutT_ZGD>> {
    using ExecPolicy =
      KernelPolicy<
        CudaKernel<
          For<2, cuda_block_x_loop, // group
            For<0, cuda_block_y_loop, // moment
              For<3, cuda_thread_x_loop, // zone
                For<1, seq_exec, // direction
                  Lambda<0>
                >
              >
            >
          >
        >
      >;
};
#endif // KRIPKE_USE_CUDA


}
}

#endif
