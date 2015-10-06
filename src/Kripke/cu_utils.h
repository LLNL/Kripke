#ifndef KRIPKE_CU_UTILS__
#define KRIPKE_CU_UTILS__

#ifdef KRIPKE_USE_CUDA

 #ifdef KRIPKE_USE_CUBLAS
 /* Using updated (v2) interfaces to cublas and cusparse */
 #include <cublas_v2.h>
 #endif



int get_cudaGetDeviceCount();
void set_cudaSetDevice(int id);
void do_cudaDeviceSynchronize();

void set_cudaMemZeroAsync( void *ptr,  size_t size);



void * get_cudaMallocManaged(size_t size);
void * get_cudaMalloc(size_t size);
void * get_cudaMallocHost(size_t size);

void  do_cudaMemcpyH2D( void *dst, void * src,  size_t size);
void  do_cudaMemcpyD2H( void *dst, void * src,  size_t size);

void  do_cudaMemcpyH2D_Async( void *dst, void * src,  size_t size);
void  do_cudaMemcpyD2H_Async( void *dst, void * src,  size_t size);

#ifdef KRIPKE_USE_CUBLAS
cublasHandle_t get_cublasHandle();
#endif

#endif

#endif
