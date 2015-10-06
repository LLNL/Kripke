
#include <stdio.h>

#include <cuda.h>
#include <cuda_runtime.h>


#include "cu_utils.h"


int get_cudaGetDeviceCount(){
   int ndev;
   cudaGetDeviceCount(&ndev);
   return  ndev;
}

void set_cudaSetDevice(int id){
   cudaSetDevice(id);
}

void do_cudaDeviceSynchronize(){
   cudaDeviceSynchronize();
} 	


void set_cudaMemZeroAsync( void *ptr,  size_t size){
   cudaMemsetAsync(ptr,0,size);
}


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

void  do_cudaMemcpyH2D_Async( void *dst, void * src,  size_t size){
  cudaMemcpyAsync(dst, src, size, cudaMemcpyHostToDevice);
}

void  do_cudaMemcpyD2H_Async( void *dst, void * src,  size_t size){
  cudaMemcpyAsync(dst, src, size, cudaMemcpyDeviceToHost);
}

#ifdef KRIPKE_USE_CUBLAS
cublasHandle_t get_cublasHandle(){
     static cublasHandle_t handle;
     static int handle_FLAG = 0;
     if (handle_FLAG==0){
        cublasCreate(&handle);
        handle_FLAG=1;
     }      
     return handle;
}
#endif


