#ifndef KRIPKE_CU_UTILS__
#define KRIPKE_CU_UTILS__

#ifdef KRIPKE_USE_CUDA

int get_cudaGetDeviceCount();
void set_cudaSetDevice(int id);

void set_cudaMemZeroAsync( void *ptr,  size_t size);



void * get_cudaMallocManaged(size_t size);
void * get_cudaMalloc(size_t size);
void * get_cudaMallocHost(size_t size);

void  do_cudaMemcpyH2D( void *dst, void * src,  size_t size);
void  do_cudaMemcpyD2H( void *dst, void * src,  size_t size);

void  do_cudaMemcpyH2D_Async( void *dst, void * src,  size_t size);
void  do_cudaMemcpyD2H_Async( void *dst, void * src,  size_t size);


#endif

#endif
