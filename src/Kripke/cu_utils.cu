
#include <stdio.h>

#include <cuda.h>
#include <cuda_runtime.h>


#include "cu_utils.h"


void * get_cudaMallocManaged(size_t size){
   void *data;
   cudaMallocManaged(&data, size);
   return data;
}



void * get_cudaMalloc(size_t size){
   void *data;
   cudaMalloc(&data, size);
   return data;
}


void * get_cudaMallocHost(size_t size){
   void *data;
   cudaMallocHost(&data, size);
   return data;
}


void  do_cudaMemcpyH2D( void *dst, void * src,  size_t size){
  cudaMemcpy(dst, src, size, cudaMemcpyHostToDevice);
}

void  do_cudaMemcpyD2H( void *dst, void * src,  size_t size){
  cudaMemcpy(dst, src, size, cudaMemcpyDeviceToHost);
}



