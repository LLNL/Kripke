#include <stdio.h>

#include <cuda.h>
#include <cuda_runtime.h>


__global__ void set_array (double *a, double value, int len)
{
    int i;
    int totalThreads = gridDim.x * blockDim.x;
    int ctaStart = blockDim.x * blockIdx.x;
    for (i = ctaStart + threadIdx.x; i < len; i += totalThreads) {
        a[i] += value;
    }
}



void set_array_test_hello(int N){
  printf("set_array_test_hello: Hello\n");

  double *h_data;
  double *d_data;
  
  h_data = new double [N];
  cudaMalloc((void **) &d_data, N*sizeof(double));

  for (int i = 0; i < N; ++i) 
    h_data[i] = i*0.1;
  
  cudaMemcpy(d_data,h_data, N*sizeof(double), cudaMemcpyHostToDevice);

  
  set_array<<<256,((N+255)/256)>>>(d_data, 2.3, N);

  cudaDeviceSynchronize();

  cudaMemcpy(h_data,d_data, N*sizeof(double), cudaMemcpyDeviceToHost);

  printf("h_data[0] = %g\n",h_data[0]);

  cudaFree(d_data);
  delete[] h_data;


}




