//
// Copyright (c) 2014-19, Lawrence Livermore National Security, LLC
// and Kripke project contributors. See the COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
//

#ifndef KRIPKE_ARCH_POPULATION
#define KRIPKE_ARCH_POPULATION

#include <Kripke.h>
#include <Kripke/VarTypes.h>

namespace Kripke {
namespace Arch {


template<typename AL>
struct Policy_Population;

template<>
struct Policy_Population<ArchLayoutT<ArchT_Sequential, LayoutT_DGZ>>{
  using ReducePolicy = seq_reduce;
  
  using ExecPolicy = 
    KernelPolicy<
      For<0, loop_exec, // direction
        For<1, loop_exec, // group
          For<2, loop_exec, // zone
            Lambda<0>
          >
        >
      >
    >;
};

template<>
struct Policy_Population<ArchLayoutT<ArchT_Sequential, LayoutT_DZG>>{
  using ReducePolicy = seq_reduce;
  
  using ExecPolicy = 
    KernelPolicy<
      For<0, loop_exec, // direction
        For<2, loop_exec, // zone
          For<1, loop_exec, // group
            Lambda<0>
          >
        >
      >
    >;
};

template<>
struct Policy_Population<ArchLayoutT<ArchT_Sequential, LayoutT_GDZ>>{
  using ReducePolicy = seq_reduce;
  
  using ExecPolicy = 
    KernelPolicy<
      For<1, loop_exec, // group
        For<0, loop_exec, // direction
          For<2, loop_exec, // zone
            Lambda<0>
          >
        >
      >
    >;
};


template<>
struct Policy_Population<ArchLayoutT<ArchT_Sequential, LayoutT_GZD>>{
  using ReducePolicy = seq_reduce;
  
  using ExecPolicy = 
    KernelPolicy<
      For<1, loop_exec, // group
        For<2, loop_exec, // zone
          For<0, loop_exec, // direction
            Lambda<0>
          >
        >
      >
    >;
};

template<>
struct Policy_Population<ArchLayoutT<ArchT_Sequential, LayoutT_ZDG>>{
  using ReducePolicy = seq_reduce;
  
  using ExecPolicy = 
    KernelPolicy<
      For<2, loop_exec, // zone
        For<0, loop_exec, // direction
          For<1, loop_exec, // group
            Lambda<0>
          >
        >
      >
    >;
};

template<>
struct Policy_Population<ArchLayoutT<ArchT_Sequential, LayoutT_ZGD>>{
  using ReducePolicy = seq_reduce;
  
  using ExecPolicy = 
    KernelPolicy<
      For<2, loop_exec, // zone
        For<1, loop_exec, // group
          For<0, loop_exec, // direction
            Lambda<0>
          >
        >
      >
    >;
};


#ifdef KRIPKE_USE_OPENMP

template<>
struct Policy_Population<ArchLayoutT<ArchT_OpenMP, LayoutT_DGZ>>{
  using ReducePolicy = omp_reduce;

  using ExecPolicy =
    KernelPolicy<
      Collapse<omp_parallel_collapse_exec, ArgList<0,1>, // Direction Group
        For<2, loop_exec, // Zone
          Lambda<0>
        >
      >
    >;
};

template<>
struct Policy_Population<ArchLayoutT<ArchT_OpenMP, LayoutT_DZG>>{
  using ReducePolicy = omp_reduce;

  using ExecPolicy =
    KernelPolicy<
      Collapse<omp_parallel_collapse_exec, ArgList<0,2>, // Direction Zone
        For<1, loop_exec, // Group
          Lambda<0>
        >
      >
    >;
};

template<>
struct Policy_Population<ArchLayoutT<ArchT_OpenMP, LayoutT_GDZ>>{
  using ReducePolicy = omp_reduce;

  using ExecPolicy =
    KernelPolicy<
      Collapse<omp_parallel_collapse_exec, ArgList<1,0>, // Group Direction
        For<2, loop_exec, // Zone
          Lambda<0>
        >
      >
    >;
};


template<>
struct Policy_Population<ArchLayoutT<ArchT_OpenMP, LayoutT_GZD>>{
  using ReducePolicy = omp_reduce;

  using ExecPolicy =
    KernelPolicy<
      Collapse<omp_parallel_collapse_exec, ArgList<1,2>, // Group Zone
        For<0, loop_exec, // Direction
          Lambda<0>
        >
      >
    >;
};

template<>
struct Policy_Population<ArchLayoutT<ArchT_OpenMP, LayoutT_ZDG>>{
  using ReducePolicy = omp_reduce;

  using ExecPolicy =
    KernelPolicy<
      Collapse<omp_parallel_collapse_exec, ArgList<2,0>, // Zone Direction
        For<1, loop_exec, // Group
          Lambda<0>
        >
      >
    >;
};

template<>
struct Policy_Population<ArchLayoutT<ArchT_OpenMP, LayoutT_ZGD>>{
  using ReducePolicy = omp_reduce;

  using ExecPolicy =
    KernelPolicy<
      Collapse<omp_parallel_collapse_exec, ArgList<2,1>, // Zone Group
        For<0, loop_exec, // Direction
          Lambda<0>
        >
      >
    >;
};


#endif // KRIPKE_USE_OPENMP



#ifdef KRIPKE_USE_CUDA
template<>
struct Policy_Population<ArchLayoutT<ArchT_CUDA, LayoutT_DGZ>>{
  using ReducePolicy = cuda_reduce;

  using ExecPolicy =
    KernelPolicy<
      CudaKernel<
        Tile<2, tile_fixed<32>, cuda_block_x_loop, // blocks of 32 zones
          For<0, cuda_thread_z_loop, // direction
            For<1, cuda_thread_y_loop, // group
              For<2, cuda_thread_x_direct, // zone
                Lambda<0>
              >
            >
          >
        >
      >
    >;
};

template<>
struct Policy_Population<ArchLayoutT<ArchT_CUDA, LayoutT_DZG>>{
  using ReducePolicy = cuda_reduce;

  using ExecPolicy =
    KernelPolicy<
      CudaKernel<
        Tile<2, tile_fixed<32>, cuda_block_x_loop, // blocks of 32 zones
          For<0, cuda_thread_z_loop, // direction
            For<2, cuda_thread_y_direct, // zone
              For<1, cuda_thread_x_loop, // group
                Lambda<0>
              >
            >
          >
        >
      >
    >;

};

template<>
struct Policy_Population<ArchLayoutT<ArchT_CUDA, LayoutT_GDZ>>{
  using ReducePolicy = cuda_reduce;

  using ExecPolicy =
    KernelPolicy<
      CudaKernel<
        Tile<2, tile_fixed<32>, cuda_block_x_loop, // blocks of 32 zones
          For<1, cuda_thread_z_loop, // group
            For<0, cuda_thread_y_loop, // direction
              For<2, cuda_thread_x_direct, // zone
                Lambda<0>
              >
            >
          >
        >
      >
    >;

};


template<>
struct Policy_Population<ArchLayoutT<ArchT_CUDA, LayoutT_GZD>>{
  using ReducePolicy = cuda_reduce;

  using ExecPolicy =
    KernelPolicy<
      CudaKernel<
        Tile<2, tile_fixed<32>, cuda_block_x_loop, // blocks of 32 zones
          For<1, cuda_thread_z_loop, // group
            For<2, cuda_thread_y_direct, // zone
              For<0, cuda_thread_x_loop, // direction
                Lambda<0>
              >
            >
          >
        >
      >
    >;

};

template<>
struct Policy_Population<ArchLayoutT<ArchT_CUDA, LayoutT_ZDG>>{
  using ReducePolicy = cuda_reduce;

  using ExecPolicy =
    KernelPolicy<
      CudaKernel<
        Tile<2, tile_fixed<32>, cuda_block_x_loop, // blocks of 32 zones
          For<2, cuda_thread_z_direct, // zone
            For<0, cuda_thread_y_loop, // direction
              For<1, cuda_thread_x_loop, // group
                Lambda<0>
              >
            >
          >
        >
      >
    >;
};

template<>
struct Policy_Population<ArchLayoutT<ArchT_CUDA, LayoutT_ZGD>>{
  using ReducePolicy = cuda_reduce;

  using ExecPolicy =
    KernelPolicy<
      CudaKernel<
        Tile<2, tile_fixed<32>, cuda_block_x_loop, // blocks of 32 zones
          For<2, cuda_thread_z_direct, // zone
            For<1, cuda_thread_y_loop, // group
              For<0, cuda_thread_x_loop, // direction
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
