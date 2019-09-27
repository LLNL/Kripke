//
// Copyright (c) 2014-19, Lawrence Livermore National Security, LLC
// and Kripke project contributors. See the COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
//

#ifndef KRIPKE_ARCH_LPLUSTIMES
#define KRIPKE_ARCH_LPLUSTIMES

#include <Kripke.h>
#include <Kripke/VarTypes.h>

namespace Kripke {
namespace Arch {


template<typename AL>
struct Policy_LPlusTimes;

template<>
struct Policy_LPlusTimes<ArchLayoutT<ArchT_Sequential, LayoutT_DGZ>> {
  using ExecPolicy =
    KernelPolicy<
      For<0, loop_exec, // Direction
        For<1, loop_exec, // Moment
          For<2, loop_exec, // Group
            For<3, loop_exec, // Zone
              Lambda<0>
            >
          >
        >
      >
    >;
};

template<>
struct Policy_LPlusTimes<ArchLayoutT<ArchT_Sequential, LayoutT_DZG>> {
  using ExecPolicy =
    KernelPolicy<
      For<0, loop_exec, // Direction
        For<1, loop_exec, // Moment
          For<3, loop_exec, // Zone
            For<2, loop_exec, // Group
              Lambda<0>
            >
          >
        >
      >
    >;
};

template<>
struct Policy_LPlusTimes<ArchLayoutT<ArchT_Sequential, LayoutT_GDZ>> {
  using ExecPolicy =
    KernelPolicy<
      For<2, loop_exec, // Group
        For<0, loop_exec, // Direction
          For<1, loop_exec, // Moment
            For<3, loop_exec, // Zone
              Lambda<0>
            >
          >
        >
      >
    >;
};


template<>
struct Policy_LPlusTimes<ArchLayoutT<ArchT_Sequential, LayoutT_GZD>> {
  using ExecPolicy =
    KernelPolicy<
      For<2, loop_exec, // Group
        For<3, loop_exec, // Zone
          For<0, loop_exec, // Direction
            For<1, loop_exec, // Moment
              Lambda<0>
            >
          >
        >
      >
    >;
};

template<>
struct Policy_LPlusTimes<ArchLayoutT<ArchT_Sequential, LayoutT_ZDG>> {
  using ExecPolicy =
    KernelPolicy<
      For<3, loop_exec, // Zone
        For<0, loop_exec, // Direction
          For<1, loop_exec, // Moment
            For<2, loop_exec, // Group
              Lambda<0>
            >
          >
        >
      >
    >;
};



template<>
struct Policy_LPlusTimes<ArchLayoutT<ArchT_Sequential, LayoutT_ZGD>> {
  using ExecPolicy =
    KernelPolicy<
      For<3, loop_exec, // Zone
        For<2, loop_exec, // Group
          For<0, loop_exec, // Direction
            For<1, loop_exec, // Moment
              Lambda<0>
            >
          >
        >
      >
    >;
};



#ifdef KRIPKE_USE_OPENMP

template<>
struct Policy_LPlusTimes<ArchLayoutT<ArchT_OpenMP, LayoutT_DGZ>> {
  using ExecPolicy =
    KernelPolicy<
      Collapse<omp_parallel_collapse_exec, ArgList<0,2>, // Direction, Group
        For<1, loop_exec, // Moment
          For<3, loop_exec, // Zone
            Lambda<0>
          >
        >
      >
    >;
};

template<>
struct Policy_LPlusTimes<ArchLayoutT<ArchT_OpenMP, LayoutT_DZG>> {
  using ExecPolicy =
    KernelPolicy<
      Collapse<omp_parallel_collapse_exec, ArgList<0,3>, // Direction, Zone
        For<1, loop_exec, // Moment
          For<2, loop_exec, // Group
            Lambda<0>
          >
        >
      >
    >;
};

template<>
struct Policy_LPlusTimes<ArchLayoutT<ArchT_OpenMP, LayoutT_GDZ>> {
  using ExecPolicy =
    KernelPolicy<
      Collapse<omp_parallel_collapse_exec, ArgList<2,0>, // Group, Direciton
        For<1, loop_exec, // Moment
          For<3, loop_exec, // Zone
            Lambda<0>
          >
        >
      >
    >;
};


template<>
struct Policy_LPlusTimes<ArchLayoutT<ArchT_OpenMP, LayoutT_GZD>> {
  using ExecPolicy =
    KernelPolicy<
      Collapse<omp_parallel_collapse_exec, ArgList<2,3,0>, // Group, Zone, Direciton
        For<1, loop_exec, // Moment
          Lambda<0>
        >
      >
    >;
};

template<>
struct Policy_LPlusTimes<ArchLayoutT<ArchT_OpenMP, LayoutT_ZDG>> {
  using ExecPolicy =
    KernelPolicy<
      Collapse<omp_parallel_collapse_exec, ArgList<3,0>, // Zone, Direction
        For<1, loop_exec, // Moment
          For<2, loop_exec, // Group
            Lambda<0>
          >
        >
      >
    >;
};



template<>
struct Policy_LPlusTimes<ArchLayoutT<ArchT_OpenMP, LayoutT_ZGD>> {
  using ExecPolicy =
    KernelPolicy<
      Collapse<omp_parallel_collapse_exec, ArgList<3,2,0>, // Zone, Group, Direction
        For<1, loop_exec, // Moment
          Lambda<0>
        >
      >
    >;
};
#endif // KRIPKE_USE_OPENMP



#ifdef KRIPKE_USE_CUDA

template<>
struct Policy_LPlusTimes<ArchLayoutT<ArchT_CUDA, LayoutT_DGZ>> {
  using ExecPolicy =
    KernelPolicy<
      CudaKernel<
        For<0, cuda_block_x_loop, // Direction
          For<2, cuda_block_y_loop, // group
            For<3, cuda_thread_x_loop, // zone
              For<1, seq_exec, // Moment
                Lambda<0>
              >
            >
          >
        >
      >
    >;
};

template<>
struct Policy_LPlusTimes<ArchLayoutT<ArchT_CUDA, LayoutT_DZG>> {
    using ExecPolicy =
      KernelPolicy<
        CudaKernel<
          For<0, cuda_block_x_loop, // Direction
            For<2, cuda_block_y_loop, // group
              For<3, cuda_thread_x_loop, // zone
                For<1, seq_exec, // Moment
                  Lambda<0>
                >
              >
            >
          >
        >
      >;
};

template<>
struct Policy_LPlusTimes<ArchLayoutT<ArchT_CUDA, LayoutT_GDZ>> {
    using ExecPolicy =
      KernelPolicy<
        CudaKernel<
          For<0, cuda_block_x_loop, // Direction
            For<2, cuda_block_y_loop, // group
              For<3, cuda_thread_x_loop, // zone
                For<1, seq_exec, // Moment
                  Lambda<0>
                >
              >
            >
          >
        >
      >;
};


template<>
struct Policy_LPlusTimes<ArchLayoutT<ArchT_CUDA, LayoutT_GZD>> {
    using ExecPolicy =
      KernelPolicy<
        CudaKernel<
          For<0, cuda_block_x_loop, // Direction
            For<2, cuda_block_y_loop, // group
              For<3, cuda_thread_x_loop, // zone
                For<1, seq_exec, // Moment
                  Lambda<0>
                >
              >
            >
          >
        >
      >;
};

template<>
struct Policy_LPlusTimes<ArchLayoutT<ArchT_CUDA, LayoutT_ZDG>> {
    using ExecPolicy =
      KernelPolicy<
        CudaKernel<
          For<0, cuda_block_x_loop, // Direction
            For<2, cuda_block_y_loop, // group
              For<3, cuda_thread_x_loop, // zone
                For<1, seq_exec, // Moment
                  Lambda<0>
                >
              >
            >
          >
        >
      >;
};



template<>
struct Policy_LPlusTimes<ArchLayoutT<ArchT_CUDA, LayoutT_ZGD>> {
    using ExecPolicy =
      KernelPolicy<
        CudaKernel<
          For<0, cuda_block_x_loop, // Direction
            For<2, cuda_block_y_loop, // group
              For<3, cuda_thread_x_loop, // zone
                For<1, seq_exec, // Moment
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
