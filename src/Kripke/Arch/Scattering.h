//
// Copyright (c) 2014-19, Lawrence Livermore National Security, LLC
// and Kripke project contributors. See the COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
//

#ifndef KRIPKE_ARCH_SCATTERING
#define KRIPKE_ARCH_SCATTERING

#include <Kripke.h>
#include <Kripke/VarTypes.h>

namespace Kripke {
namespace Arch {

template<typename AL>
struct Policy_Scattering;

template<>
struct Policy_Scattering<ArchLayoutT<ArchT_Sequential, LayoutT_DGZ>> {
  using ExecPolicy =
      KernelPolicy<
        For<0, loop_exec, // moment
          For<1, loop_exec, // dst group
            For<2, loop_exec, // src group
              For<3, loop_exec, // zone
                Lambda<0>
              >
            >
          >
        >
      >;
};

template<>
struct Policy_Scattering<ArchLayoutT<ArchT_Sequential, LayoutT_DZG>> {
  using ExecPolicy =
      KernelPolicy<
        For<0, loop_exec, // moment
          For<3, loop_exec, // zone
            For<1, loop_exec, // dst group
              For<2, loop_exec, // src group
                Lambda<0>
              >
            >
          >
        >
      >;
};


template<>
struct Policy_Scattering<ArchLayoutT<ArchT_Sequential, LayoutT_GDZ>> {
  using ExecPolicy =
      KernelPolicy<
        For<1, loop_exec, // dst group
          For<2, loop_exec, // src group
            For<0, loop_exec, // moment
              For<3, loop_exec, // zone
                Lambda<0>
              >
            >
          >
        >
      >;
};


template<>
struct Policy_Scattering<ArchLayoutT<ArchT_Sequential, LayoutT_GZD>> {
  using ExecPolicy =
      KernelPolicy<
        For<1, loop_exec, // dst group
          For<2, loop_exec, // src group
            For<3, loop_exec, // zone
              For<0, loop_exec, // moment
                Lambda<0>
              >
            >
          >
        >
      >;
};


template<>
struct Policy_Scattering<ArchLayoutT<ArchT_Sequential, LayoutT_ZDG>> {
  using ExecPolicy =
      KernelPolicy<
        For<3, loop_exec, // zone
          For<0, loop_exec, // moment
            For<1, loop_exec, // dst group
              For<2, loop_exec, // src group
                Lambda<0>
              >
            >
          >
        >
      >;
};


template<>
struct Policy_Scattering<ArchLayoutT<ArchT_Sequential, LayoutT_ZGD>> {
  using ExecPolicy =
      KernelPolicy<
        For<3, loop_exec, // zone
          For<1, loop_exec, // dst group
            For<2, loop_exec, // src group
              For<0, loop_exec, // moment
                Lambda<0>
              >
            >
          >
        >
      >;
};



#ifdef KRIPKE_USE_OPENMP
template<>
struct Policy_Scattering<ArchLayoutT<ArchT_OpenMP, LayoutT_DGZ>> {
  using ExecPolicy =
    KernelPolicy<
      Collapse<omp_parallel_collapse_exec, ArgList<0,1>, // Moment, DstGrp
        For<2, loop_exec, // SrcGrp
          For<3, loop_exec, // Zone
            Lambda<0>
          >
        >
      >
    >;
};

template<>
struct Policy_Scattering<ArchLayoutT<ArchT_OpenMP, LayoutT_DZG>> {
  using ExecPolicy =
    KernelPolicy<
      Collapse<omp_parallel_collapse_exec, ArgList<0,3,1>, // Moment, Zone, DstGrp
        For<2, loop_exec, // SrcGrp
          Lambda<0>
        >
      >
    >;
};


template<>
struct Policy_Scattering<ArchLayoutT<ArchT_OpenMP, LayoutT_GDZ>> {
  using ExecPolicy =
    KernelPolicy<
      Collapse<omp_parallel_collapse_exec, ArgList<1,0>, // DstGrp, Moment
        For<2, loop_exec, // SrcGrp
          For<3, loop_exec, // Zone
            Lambda<0>
          >
        >
      >
    >;
};


template<>
struct Policy_Scattering<ArchLayoutT<ArchT_OpenMP, LayoutT_GZD>> {
  using ExecPolicy =
    KernelPolicy<
      Collapse<omp_parallel_collapse_exec, ArgList<1,3>, // DstGrp, Zone
        For<2, loop_exec, // SrcGrp
          For<0, loop_exec, // Moment
            Lambda<0>
          >
        >
      >
    >;
};


template<>
struct Policy_Scattering<ArchLayoutT<ArchT_OpenMP, LayoutT_ZDG>> {
  using ExecPolicy =
    KernelPolicy<
      Collapse<omp_parallel_collapse_exec, ArgList<3,0,1>, // Zone, Moment, DstGrp
        For<2, loop_exec, // SrcGrp
          Lambda<0>
        >
      >
    >;
};


template<>
struct Policy_Scattering<ArchLayoutT<ArchT_OpenMP, LayoutT_ZGD>> {
  using ExecPolicy =
    KernelPolicy<
      Collapse<omp_parallel_collapse_exec, ArgList<3,1>, // Zone, DstGrp
        For<2, loop_exec, // SrcGrp
          For<0, loop_exec, // Moment
            Lambda<0>
          >
        >
      >
    >;
};
#endif // KRIPKE_USE_OPENMP


#ifdef KRIPKE_USE_CUDA
template<>
struct Policy_Scattering<ArchLayoutT<ArchT_CUDA, LayoutT_DGZ>> {
  using ExecPolicy =
    KernelPolicy<
      CudaKernel<
        For<0, cuda_block_x_loop, // moment
          For<1, cuda_block_y_loop, // DstGrp
            For<3, cuda_thread_x_loop, // zone
              For<2, seq_exec, // SrcGrp
                Lambda<0>
              >
            >
          >
        >
      >
    >;
};

template<>
struct Policy_Scattering<ArchLayoutT<ArchT_CUDA, LayoutT_DZG>> {
    using ExecPolicy =
      KernelPolicy<
        CudaKernel<
          For<0, cuda_block_x_loop, // moment
            For<1, cuda_block_y_loop, // DstGrp
              For<3, cuda_thread_x_loop, // zone
                For<2, seq_exec, // SrcGrp
                  Lambda<0>
                >
              >
            >
          >
        >
      >;
};


template<>
struct Policy_Scattering<ArchLayoutT<ArchT_CUDA, LayoutT_GDZ>> {
    using ExecPolicy =
      KernelPolicy<
        CudaKernel<
          For<0, cuda_block_x_loop, // moment
            For<1, cuda_block_y_loop, // DstGrp
              For<3, cuda_thread_x_loop, // zone
                For<2, seq_exec, // SrcGrp
                  Lambda<0>
                >
              >
            >
          >
        >
      >;
};


template<>
struct Policy_Scattering<ArchLayoutT<ArchT_CUDA, LayoutT_GZD>> {
    using ExecPolicy =
      KernelPolicy<
        CudaKernel<
          For<0, cuda_block_x_loop, // moment
            For<1, cuda_block_y_loop, // DstGrp
              For<3, cuda_thread_x_loop, // zone
                For<2, seq_exec, // SrcGrp
                  Lambda<0>
                >
              >
            >
          >
        >
      >;
};


template<>
struct Policy_Scattering<ArchLayoutT<ArchT_CUDA, LayoutT_ZDG>> {
    using ExecPolicy =
      KernelPolicy<
        CudaKernel<
          For<0, cuda_block_x_loop, // moment
            For<1, cuda_block_y_loop, // DstGrp
              For<3, cuda_thread_x_loop, // zone
                For<2, seq_exec, // SrcGrp
                  Lambda<0>
                >
              >
            >
          >
        >
      >;
};


template<>
struct Policy_Scattering<ArchLayoutT<ArchT_CUDA, LayoutT_ZGD>> {
    using ExecPolicy =
      KernelPolicy<
        CudaKernel<
          For<0, cuda_block_x_loop, // moment
            For<1, cuda_block_y_loop, // DstGrp
              For<3, cuda_thread_x_loop, // zone
                For<2, seq_exec, // SrcGrp
                  Lambda<0>
                >
              >
            >
          >
        >
      >;
};
#endif //KRIPKE_USE_CUDA


}
}

#endif
