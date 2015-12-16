#ifndef RAJA_CUDA_COMMON_H
#define RAJA_CUDA_COMMON_H

#ifdef RAJA_USE_CUDA


struct CudaDim {
  dim3 num_threads;
  dim3 num_blocks;
  
  __host__ __device__ void print(void){
    printf("<<< (%d,%d,%d), (%d,%d,%d) >>>\n",
      num_blocks.x, num_blocks.y, num_blocks.z,
      num_threads.x, num_threads.y, num_threads.z);
  }
};


struct Dim3x {
  __host__ __device__ inline unsigned int &operator()(dim3 &dim){
    return dim.x;
  }
  
  __host__ __device__ inline unsigned int operator()(dim3 const &dim){
    return dim.x;
  }
};


struct Dim3y {
  __host__ __device__ inline unsigned int &operator()(dim3 &dim){
    return dim.y;
  }
  
  __host__ __device__ inline unsigned int operator()(dim3 const &dim){
    return dim.y;
  }
};

struct Dim3z {
  __host__ __device__ inline unsigned int &operator()(dim3 &dim){
    return dim.z;
  }
  
  __host__ __device__ inline unsigned int operator()(dim3 const &dim){
    return dim.z;
  }
};

template<typename POL>
struct CudaPolicy {};

template<typename VIEWDIM, int threads_per_block>
struct CudaThreadBlock {
  int begin;
  int end;
  
  VIEWDIM view;
  
  CudaThreadBlock(int begin0, int end0) : begin(begin0), end(end0){}

  __device__ inline int operator()(void){
    
    int idx = begin + view(blockIdx) * threads_per_block + view(threadIdx);
    if(idx >= end){
      idx = -1;
    }
    return idx;
  }
  
  void inline setDims(CudaDim &dims){
    int n = end-begin;
    if(n < threads_per_block){
      view(dims.num_threads) = n;
      view(dims.num_blocks) = 1;
    }
    else{
      view(dims.num_threads) = threads_per_block;
      
      int blocks = n / threads_per_block;
      if(n % threads_per_block){
        ++ blocks;
      }
      view(dims.num_blocks) = blocks;
    }
  }  
};

#endif // RAJA_USE_CUDA

#endif

