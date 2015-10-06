#include <stdio.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "Directions.h"


#ifdef KRIPKE_USE_CUBLAS
/* Using updated (v2) interfaces to cublas and cusparse */
#include <cublas_v2.h>
#include "cu_utils.h"
#endif



#define KRESTRICT __restrict__

#define USE_PSI_HOST_MEM

//#define USE_IJK_PLANE_HOST_MEM

#define CU_TIMING__

#define MAX(a,b) ((a<b)?b:a)
#define MIN(a,b) ((a>b)?b:a)

//Macro for checking cuda errors following a cuda launch or api call
#define cudaCheckError() {                                              \
  cudaError_t e=cudaGetLastError();                                     \
  if(e!=cudaSuccess) {                                                  \
    printf("Cuda failure %s:%d: '%s'\n",__FILE__,__LINE__,cudaGetErrorString(e)); \
    exit(EXIT_FAILURE);                                                 \
  }                                                                     \
}


int cuda_LTimes_GZD(double *phi, double *psi, double *ell,
                    int num_groups_zones, int num_local_directions, int nidx, int group0);


__global__ void LTimes_GZD(double *phi, double *psi, double *ell,
                    int num_groups_zones, int num_local_directions, int nidx, int group0);


int cuda_LTimes_GZD(double *phi, double *psi, double *ell,
                    int num_groups_zones, int num_local_directions, int nidx, int group0){

#if 1
  //call cuBLAS
  double ONE = 1.0;
  double ZERO = 0.0;
  cublasOperation_t transa = CUBLAS_OP_N;
  cublasOperation_t transb = CUBLAS_OP_N;

  cublasHandle_t handle;
  handle = get_cublasHandle();


  cudaCheckError();

  cublasDgemm(handle,
              transa, transb,
              nidx, num_groups_zones, num_local_directions,
              &ONE,
              (const double *) ell,  nidx,
              (const double  *) psi, num_local_directions,
              &ONE,
              phi, nidx);

  

  cudaCheckError();
#else
  LTimes_GZD<<<num_groups_zones,32>>>(phi,psi,ell,num_groups_zones,num_local_directions,nidx,group0);

#endif

  cudaDeviceSynchronize();
  return 0;

}

__global__ void LTimes_GZD(double *phi, double *psi, double *ell,
                    int num_groups_zones, int num_local_directions, int nidx, int group0){

        
        for (int i = blockIdx.x; i < num_groups_zones; i+=gridDim.x)
          for (int j = 0; j < num_local_directions; ++j)
            for (int k = threadIdx.x; k < nidx; k += blockDim.x)
              phi[i*nidx + k] += ell[j*nidx + k] * psi[i*num_local_directions + j];

}


int cuda_sweep_GZD( double *d_rhs, double *h_phi, double *h_psi, double *d_sigt, Directions *d_direction,
                    double *h_i_plane, double *h_j_plane, double *h_k_plane,
                    int *d_ii_jj_kk_z_idx, int *h_offset,  int *d_offset, double *d_dx, double *d_dy, double *d_dz,
                    int num_zones, int num_directions, int num_groups,
                    int local_imax, int local_jmax, int local_kmax, int Nslices,
                    int start_i, int start_j, int start_k, 
                    int end_i, int end_j, int end_k,
                    int inc_i, int inc_j, int inc_k);


__global__  void sweep_ijk_GDZ(const int num_zones,  const int num_directions, const int num_groups,
                               const int local_imax, const int local_jmax,     const int local_kmax,
                               const double * KRESTRICT d_dx, 
                               const double * KRESTRICT d_dy, 
                               const double * KRESTRICT d_dz, 
                               double * KRESTRICT d_rhs, 
                               double * KRESTRICT d_psi, 
                               const double * KRESTRICT d_sigt, 
                               const Directions * KRESTRICT d_direction,
                               double * KRESTRICT d_i_plane, 
                               double * KRESTRICT d_j_plane, 
                               double * KRESTRICT d_k_plane,
                               const int start_i, const int start_j, const int start_k, 
                               const int end_i, const int end_j, const int end_k,
                               const int inc_i, const int inc_j, const int inc_k);


#define I_PLANE_INDEX(j, k) ((k)*(local_jmax) + (j))
#define J_PLANE_INDEX(i, k) ((k)*(local_imax) + (i))
#define K_PLANE_INDEX(i, j) ((j)*(local_imax) + (i))
#define Zonal_INDEX(i, j, k) ((i) + (local_imax)*(j) \
  + (local_imax)*(local_jmax)*(k))


int cuda_sweep_GZD( double *d_rhs, double *h_phi, double *h_psi, double *d_sigt, Directions *d_direction,
                    double *h_i_plane, double *h_j_plane, double *h_k_plane,
                    int *d_ii_jj_kk_z_idx, int *h_offset,  int *d_offset, double *d_dx, double *d_dy, double *d_dz,
                    int num_zones, int num_directions, int num_groups,
                    int local_imax, int local_jmax, int local_kmax, int Nslices,
                    int start_i, int start_j, int start_k, 
                    int end_i, int end_j, int end_k,
                    int inc_i, int inc_j, int inc_k){

    dim3 numBlocks(num_groups);
    dim3 threadsPerBlock(32);

    double *d_phi = h_phi;  //use zero_copy

  #ifdef CU_TIMING
  float time_ms, time_s;
  cudaEvent_t start,stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  #endif

#ifndef USE_PSI_HOST_MEM
  double *d_psi;
  int N = num_zones*num_directions*num_groups;
  #ifdef CU_TIMING
  cudaEventRecord(start);
  #endif
  cudaMalloc((void **) &d_psi, N*sizeof(double));
//  cudaMemcpy(d_psi, h_psi,   N*sizeof(double), cudaMemcpyHostToDevice);
  cudaCheckError();

  #ifdef CU_TIMING
  cudaEventRecord(stop);
  cudaDeviceSynchronize();
  cudaCheckError();
  cudaEventElapsedTime(&time_ms,start,stop);
  time_s=time_ms*.001;
  printf("GZD: time to allocated memory for PSI on device: %g [s]\n",time_s);
  #endif

#else
    double *d_psi = h_psi;  //use zero copy
#endif


    int groups_dirs = num_groups*num_directions;
 
    int i_plane_zones = local_jmax * local_kmax * groups_dirs;
    int j_plane_zones = local_imax * local_kmax * groups_dirs;
    int k_plane_zones = local_imax * local_jmax * groups_dirs;

#ifdef USE_IJK_PLANE_HOST_MEM

  double *d_i_plane = h_i_plane;
  double *d_j_plane = h_j_plane;
  double *d_k_plane = h_k_plane;

#else

  double *d_i_plane,  *d_j_plane, *d_k_plane;

  #ifdef CU_TIMING
  cudaEventRecord(start);
  #endif

  cudaMalloc((void **) &d_i_plane, i_plane_zones * sizeof(double));
  cudaMalloc((void **) &d_j_plane, j_plane_zones * sizeof(double));
  cudaMalloc((void **) &d_k_plane, k_plane_zones * sizeof(double));
  cudaMemcpy(d_i_plane, h_i_plane, i_plane_zones * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(d_j_plane, h_j_plane, j_plane_zones * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(d_k_plane, h_k_plane, k_plane_zones * sizeof(double), cudaMemcpyHostToDevice);
  cudaCheckError();

  #ifdef CU_TIMING
  cudaEventRecord(stop);
  cudaDeviceSynchronize();
  cudaCheckError();
  cudaEventElapsedTime(&time_ms,start,stop);
  time_s=time_ms*.001;
  printf("ZGD: time to copy ijk_PLANE H2D: %g [s]\n",time_s);
  #endif


#endif



    #ifdef CU_TIMING
    cudaEventRecord(start);
    #endif

    int shared_memory_size = (num_directions*12 + 12)*sizeof(double);
    sweep_ijk_GDZ<<<numBlocks,threadsPerBlock,shared_memory_size>>>(num_zones,num_directions,num_groups,local_imax,local_jmax,local_kmax,
                                                 d_dx, d_dy, d_dz, d_rhs,  d_psi, d_sigt, d_direction, 
                                                 d_i_plane, d_j_plane, d_k_plane, start_i, start_j, start_k,
                                                 end_i, end_j, end_k, inc_i, inc_j, inc_k);

    cudaCheckError();

    #ifdef CU_TIMING
    cudaEventRecord(stop);
    cudaDeviceSynchronize();
    cudaCheckError();
    cudaEventElapsedTime(&time_ms,start,stop);
    time_s=time_ms*.001;
    printf("sweep_ijk_GDZ  time: %g [s]\n",time_s);
    #endif


#ifndef USE_IJK_PLANE_HOST_MEM
    #ifdef CU_TIMING
    cudaEventRecord(start);
    #endif

    cudaMemcpy(h_i_plane, d_i_plane, i_plane_zones * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_j_plane, d_j_plane, j_plane_zones * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_k_plane, d_k_plane, k_plane_zones * sizeof(double), cudaMemcpyDeviceToHost);
    cudaFree(d_i_plane);
    cudaFree(d_j_plane);
    cudaFree(d_k_plane);

    #ifdef CU_TIMING
    cudaEventRecord(stop);
    cudaDeviceSynchronize();
    cudaCheckError();
    cudaEventElapsedTime(&time_ms,start,stop);
    time_s=time_ms*.001;
    printf("ZGD: time to copy ijk_PLANE D2H: %g [s]\n",time_s);
    #endif
#endif

#ifndef USE_PSI_HOST_MEM

  #ifdef CU_TIMING
  cudaEventRecord(start);
  #endif

  cudaMemcpy(h_psi,     d_psi, N*sizeof(double),                   cudaMemcpyDeviceToHost);
  cudaFree(d_psi);

  #ifdef CU_TIMING
  cudaEventRecord(stop);
  cudaDeviceSynchronize();
  cudaCheckError();
  cudaEventElapsedTime(&time_ms,start,stop);
  time_s=time_ms*.001;
  printf("GZD: time to copy PSI D2H: %g [s]\n",time_s);
  #endif

#endif

   return 0;
}


__global__ void sweep_ijk_GDZ( const int num_zones,  const int num_local_directions, const int num_local_groups,
                               const int local_imax, const int local_jmax,     const int local_kmax,
                               const double * KRESTRICT dx,
                               const double * KRESTRICT dy,
                               const double * KRESTRICT dz,
                               double * KRESTRICT rhs,
                               double * KRESTRICT psi,
                               const double * KRESTRICT sigt,
                               const Directions * KRESTRICT direction,
                               double * KRESTRICT i_plane,
                               double * KRESTRICT j_plane,
                               double * KRESTRICT k_plane,
                               const int start_i, const int start_j, const int start_k,
                               const int end_i, const int end_j, const int end_k,
                               const int inc_i, const int inc_j, const int inc_k){

    extern __shared__ double ss_psi_fr_g_z[];
    int i_BS = 12;
    double * __restrict__ ss_two_inv_dxi = &ss_psi_fr_g_z[i_BS*num_local_directions];
    int ii;

/*
    for (int k = start_k; k != end_k; k += inc_k) {
      const double dzk = dz[k + 1];
      const double two_dz = 2.0 / dzk;
    }

     for (int j = start_j; j != end_j; j += inc_j) {
       const double dyj = dy[j + 1];
       const double two_dy = 2.0 / dyj;      
    }


     for (int i = start_i; j != end_i; j += inc_i) {
       const double dxi = dx[i + 1];
       const double two_dx = 2.0 / dxi;
    }

*/



    for (int group = blockIdx.x; group < num_local_groups; group += gridDim.x){

      int i_plane_zones_group = local_jmax * local_kmax * group;// * groups_dirs;
      int j_plane_zones_group = local_imax * local_kmax * group;// * groups_dirs;
      int k_plane_zones_group = local_imax * local_jmax * group;// * groups_dirs;
      const double * KRESTRICT block_sigt = &sigt[group * num_zones];



/*
  offset = g*nzones*directions +
           z*directions +
           d
*/
      int group_offset = group*num_zones*num_local_directions;

    /*  Perform transport sweep of the grid 1 cell at a time.   */
      for (int k = start_k; k != end_k; k += inc_k) {
        //const double dzk = dz[k + 1];
        //const double two_dz = 2.0 / dzk;

        const double two_dz = dz[k + 1];

        //block in "i-" direction to limit use of shared memory  
        int iblock_range = (end_i-start_i)*inc_i; 
        for (int block_i = 0; block_i < (iblock_range+i_BS-1)/i_BS; block_i++){
 
 
          ii=0;
 //         for (int i = start_i; i != end_i; i += inc_i, ii++) {
          for (int iii = 0; iii < MIN(i_BS, iblock_range-block_i*i_BS); iii++, ii+=num_local_directions) {
            int i = start_i + (block_i*i_BS +iii)*inc_i;
            double * KRESTRICT psi_fr_g_z = &j_plane[(j_plane_zones_group+J_PLANE_INDEX(i, k))*num_local_directions];
            for (int d = threadIdx.x; d < num_local_directions; d += blockDim.x)
              ss_psi_fr_g_z[ii + d] = psi_fr_g_z[d];
          } 
          for (int iii = threadIdx.x; iii < MIN(i_BS, iblock_range-block_i*i_BS); iii+=blockDim.x){
              int i = start_i + (block_i*i_BS +iii)*inc_i;
              //ss_two_inv_dxi[iii] = 2.0/ dx[i + 1];
              ss_two_inv_dxi[iii] = dx[i + 1];  
          }
          __syncthreads();

          for (int j = start_j; j != end_j; j += inc_j) {
            //const double dyj = dy[j + 1];
            //const double two_dy = 2.0 / dyj;
            const double two_dy = dy[j + 1];
            double * KRESTRICT psi_lf_g_z = &i_plane[(i_plane_zones_group+I_PLANE_INDEX(j, k))*num_local_directions];
       
            //assume blockDim.x >= num_local_directions
            int i = start_i + block_i*i_BS*inc_i;
            double temp_psi_lf_g_z = psi_lf_g_z[threadIdx.x];//keep flux in registers
            
            int z = Zonal_INDEX(i, j, k);
            double * KRESTRICT psi_fr_g_z = &ss_psi_fr_g_z[0];
            double * KRESTRICT psi_bo_g_z = &k_plane[(k_plane_zones_group+ K_PLANE_INDEX(i, j) )*num_local_directions];

            double * KRESTRICT psi_g_z = &psi[group_offset + z*num_local_directions];
            double * KRESTRICT rhs_g_z = &rhs[group_offset + z*num_local_directions];


            //ii=0;
//            for (int i = start_i; i != end_i; i += inc_i, ii++) {
            for (int iii = 0; 
                 iii < MIN(i_BS, iblock_range-block_i*i_BS); 
                 iii++,z+=inc_i,psi_g_z+=(num_local_directions*inc_i), rhs_g_z+=(num_local_directions*inc_i),
                 psi_fr_g_z += num_local_directions,
                 psi_bo_g_z += (num_local_directions*inc_i) ) {
           
              double two_dx = ss_two_inv_dxi[iii];

              const  double block_sigt_z =  block_sigt[z];
              for (int d = threadIdx.x; d < num_local_directions; d += blockDim.x) {
                const double xcos = direction[d].xcos;
                const double ycos = direction[d].ycos;
                const double zcos = direction[d].zcos;

                const double xcos_dxi = xcos * two_dx;
                const double ycos_dyj = ycos * two_dy;
                const double zcos_dzk = zcos * two_dz;

                /* Calculate new zonal flux */
                double temp = 1.0/(  xcos_dxi + ycos_dyj + zcos_dzk + block_sigt_z);

                double psi_g_z_d = (rhs_g_z[d] + 
                                  temp_psi_lf_g_z  * xcos_dxi +
                                  psi_fr_g_z[d]    * ycos_dyj + 
                                  psi_bo_g_z[d]    * zcos_dzk)*temp;

                psi_g_z[d] = psi_g_z_d;

                /* Apply diamond-difference relationships */
                temp_psi_lf_g_z = 2.0 * psi_g_z_d - temp_psi_lf_g_z; //flux is kept in registers
                psi_fr_g_z[d] = 2.0 * psi_g_z_d - psi_fr_g_z[d];     //flux is kept in shared memory
                psi_bo_g_z[d] = 2.0 * psi_g_z_d - psi_bo_g_z[d];     //flux is kept in global memory 
              }
            } 
            psi_lf_g_z[threadIdx.x] = temp_psi_lf_g_z;
          }
        
          ii = 0;
//          for (int i = start_i; i != end_i; i += inc_i, ii++) {
          for (int iii = 0; iii < MIN(i_BS, iblock_range-block_i*i_BS); iii++, ii += num_local_directions) {
            int i = start_i + (block_i*i_BS + iii) * inc_i;
            double * KRESTRICT psi_fr_g_z = &j_plane[(j_plane_zones_group+J_PLANE_INDEX(i, k))*num_local_directions];
            for (int d = threadIdx.x; d < num_local_directions; d += blockDim.x)
              psi_fr_g_z[d] = ss_psi_fr_g_z[ii + d];
          }
        }// end of "for (iii=0..."
      } //end of "for (k=..."
   } // end of for (group=0;...
}
