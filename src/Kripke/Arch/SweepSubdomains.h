//
// Copyright (c) 2014-19, Lawrence Livermore National Security, LLC
// and Kripke project contributors. See the COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
//

#ifndef KRIPKE_ARCH_SWEEPSUBDOMAINS
#define KRIPKE_ARCH_SWEEPSUBDOMAINS

#include <Kripke.h>
#include <Kripke/VarTypes.h>

namespace Kripke {
namespace Arch {

template<typename AL>
struct Policy_SweepSubdomains;

template<>
struct Policy_SweepSubdomains<ArchLayoutT<ArchT_Sequential, LayoutT_DGZ>> {
  using ExecPolicy =
    KernelPolicy<
      For<0, loop_exec,  // direction
        For<1, loop_exec, // group
          For<2, loop_exec, // k
            For<3, loop_exec, // j
              For<4, loop_exec, // i
                Lambda<0>
              >
            >
          >
        >
      >
    >;
};


template<>
struct Policy_SweepSubdomains<ArchLayoutT<ArchT_Sequential, LayoutT_DZG>> {
  using ExecPolicy =
    KernelPolicy<
      For<0, loop_exec,  // direction
        For<2, loop_exec, // k
          For<3, loop_exec, // j
            For<4, loop_exec, // i
              For<1, loop_exec, // group
                Lambda<0>
              >
            >
          >
        >
      >
    >;
};


template<>
struct Policy_SweepSubdomains<ArchLayoutT<ArchT_Sequential, LayoutT_GDZ>> {
  using ExecPolicy =
    KernelPolicy<
      For<1, loop_exec, // group
        For<0, loop_exec,  // direction
          For<2, loop_exec, // k
            For<3, loop_exec, // j
              For<4, loop_exec, // i
                Lambda<0>
              >
            >
          >
        >
      >
    >;
};


template<>
struct Policy_SweepSubdomains<ArchLayoutT<ArchT_Sequential, LayoutT_GZD>> {
  using ExecPolicy =
    KernelPolicy<
      For<1, loop_exec, // group
        For<2, loop_exec, // k
          For<3, loop_exec, // j
            For<4, loop_exec, // i
              For<0, loop_exec,  // direction
                Lambda<0>
              >
            >
          >
        >
      >
    >;
};


template<>
struct Policy_SweepSubdomains<ArchLayoutT<ArchT_Sequential, LayoutT_ZDG>> {
  using ExecPolicy =
    KernelPolicy<
      For<2, loop_exec, // k
        For<3, loop_exec, // j
          For<4, loop_exec, // i
            For<0, loop_exec,  // direction
              For<1, loop_exec, // group
                Lambda<0>
              >
            >
          >
        >
      >
    >;
};


template<>
struct Policy_SweepSubdomains<ArchLayoutT<ArchT_Sequential, LayoutT_ZGD>> {
  using ExecPolicy =
    KernelPolicy<
      For<2, loop_exec, // k
        For<3, loop_exec, // j
          For<4, loop_exec, // i
            For<1, loop_exec, // group
              For<0, loop_exec,  // direction
                Lambda<0>
              >
            >
          >
        >
      >
    >;
};




#ifdef KRIPKE_USE_OPENMP


template<>
struct Policy_SweepSubdomains<ArchLayoutT<ArchT_OpenMP, LayoutT_DGZ>> {


  using ExecPolicy =
    KernelPolicy<
      Collapse<omp_parallel_collapse_exec, ArgList<0,1>, // direction, group
        For<2, loop_exec, // k
          For<3, loop_exec, // j
            For<4, loop_exec, // i
              Lambda<0>
            >
          >
        >
      >
    >;
};


template<>
struct Policy_SweepSubdomains<ArchLayoutT<ArchT_OpenMP, LayoutT_DZG>> {
  using ExecPolicy =
    KernelPolicy<
      For<0, omp_parallel_for_exec,  // direction
        For<2, loop_exec, // k
          For<3, loop_exec, // j
            For<4, loop_exec, // i
              For<1, loop_exec, // group
                Lambda<0>
              >
            >
          >
        >
      >
    >;
};


template<>
struct Policy_SweepSubdomains<ArchLayoutT<ArchT_OpenMP, LayoutT_GDZ>> {
  using ExecPolicy =
    KernelPolicy<
      Collapse<omp_parallel_collapse_exec, ArgList<1,0>, // group, direction
        For<2, loop_exec, // k
          For<3, loop_exec, // j
            For<4, loop_exec, // i
              Lambda<0>
            >
          >
        >
      >
    >;
};


template<>
struct Policy_SweepSubdomains<ArchLayoutT<ArchT_OpenMP, LayoutT_GZD>> {
  using ExecPolicy =
    KernelPolicy<
      For<1, omp_parallel_for_exec, // group
        For<2, loop_exec, // k
          For<3, loop_exec, // j
            For<4, loop_exec, // i
              For<0, loop_exec,  // direction
                Lambda<0>
              >
            >
          >
        >
      >
    >;
};


template<>
struct Policy_SweepSubdomains<ArchLayoutT<ArchT_OpenMP, LayoutT_ZDG>> {
  using ExecPolicy =
    KernelPolicy<
      Hyperplane<2, seq_exec, ArgList<3,4>, omp_parallel_collapse_exec,
        For<0, loop_exec,  // direction
          For<1, loop_exec, // group
            Lambda<0>
          >
        >
      >
    >;
};


template<>
struct Policy_SweepSubdomains<ArchLayoutT<ArchT_OpenMP, LayoutT_ZGD>> {
  using ExecPolicy =
    KernelPolicy<
      Hyperplane<2, seq_exec, ArgList<3,4>, omp_parallel_collapse_exec,
        For<1, loop_exec, // group
          For<0, loop_exec,  // direction
            Lambda<0>
          >
        >
      >
    >;
};

#endif // KRIPKE_USE_OPENMP


#ifdef KRIPKE_USE_CUDA
template<>
struct Policy_SweepSubdomains<ArchLayoutT<ArchT_CUDA, LayoutT_DGZ>> {
  using ExecPolicy =
          KernelPolicy<
            CudaKernel<
              For<0, cuda_block_x_loop,
                For<1, cuda_block_y_loop,

                  For<3, cuda_thread_y_loop,
                    For<4, cuda_thread_x_loop,
                      Hyperplane<
                        2, seq_exec, ArgList<3, 4>,

                        Lambda<0>,
                        CudaSyncThreads
                      >
                    >
                  >

                >
              >
            >
          >;
};


template<>
struct Policy_SweepSubdomains<ArchLayoutT<ArchT_CUDA, LayoutT_DZG>> {
    using ExecPolicy =
            KernelPolicy<
              CudaKernel<
                For<0, cuda_block_x_loop,
                  For<1, cuda_block_y_loop,

                    For<3, cuda_thread_y_loop,
                      For<4, cuda_thread_x_loop,
                        Hyperplane<
                          2, seq_exec, ArgList<3, 4>,

                          Lambda<0>,
                          CudaSyncThreads
                        >
                      >
                    >

                  >
                >
              >
            >;
};


template<>
struct Policy_SweepSubdomains<ArchLayoutT<ArchT_CUDA, LayoutT_GDZ>> {
    using ExecPolicy =
            KernelPolicy<
              CudaKernel<
                For<0, cuda_block_x_loop,
                  For<1, cuda_block_y_loop,

                    For<3, cuda_thread_y_loop,
                      For<4, cuda_thread_x_loop,
                        Hyperplane<
                          2, seq_exec, ArgList<3, 4>,

                          Lambda<0>,
                          CudaSyncThreads
                        >
                      >
                    >

                  >
                >
              >
            >;
};


template<>
struct Policy_SweepSubdomains<ArchLayoutT<ArchT_CUDA, LayoutT_GZD>> {
    using ExecPolicy =
            KernelPolicy<
              CudaKernel<
                For<0, cuda_block_x_loop,
                  For<1, cuda_block_y_loop,

                    For<3, cuda_thread_y_loop,
                      For<4, cuda_thread_x_loop,
                        Hyperplane<
                          2, seq_exec, ArgList<3, 4>,

                          Lambda<0>,
                          CudaSyncThreads
                        >
                      >
                    >

                  >
                >
              >
            >;
};


template<>
struct Policy_SweepSubdomains<ArchLayoutT<ArchT_CUDA, LayoutT_ZDG>> {
    using ExecPolicy =
            KernelPolicy<
              CudaKernel<
                For<0, cuda_block_x_loop,
                  For<1, cuda_block_y_loop,

                    For<3, cuda_thread_y_loop,
                      For<4, cuda_thread_x_loop,
                        Hyperplane<
                          2, seq_exec, ArgList<3, 4>,

                          Lambda<0>,
                          CudaSyncThreads
                        >
                      >
                    >

                  >
                >
              >
            >;
};


template<>
struct Policy_SweepSubdomains<ArchLayoutT<ArchT_CUDA, LayoutT_ZGD>> {
    using ExecPolicy =
            KernelPolicy<
              CudaKernel<
                For<0, cuda_block_x_loop,
                  For<1, cuda_block_y_loop,

                    For<3, cuda_thread_y_loop,
                      For<4, cuda_thread_x_loop,
                        Hyperplane<
                          2, seq_exec, ArgList<3, 4>,

                          Lambda<0>,
                          CudaSyncThreads
                        >
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
