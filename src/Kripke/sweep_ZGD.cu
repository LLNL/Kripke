
#include <stdio.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "Directions.h"


#define KRESTRICT __restrict__

#define USE_PSI_HOST_MEM

//#define CU_TIMING

#define MAX ((a<b)?b:a)
#define MIN ((a>b)?b:a)


//Macro for checking cuda errors following a cuda launch or api call
#define cudaCheckError() {                                              \
  cudaError_t e=cudaGetLastError();                                     \
  if(e!=cudaSuccess) {                                                  \
    printf("Cuda failure %s:%d: '%s'\n",__FILE__,__LINE__,cudaGetErrorString(e)); \
    exit(EXIT_FAILURE);                                                 \
  }                                                                     \
}


__global__ void  LTimes_ZGD(double *phi, double * __restrict__ psi, double * __restrict__ ell,
                                int num_zones, int num_groups, int num_local_directions, int num_local_groups, int nidx);


__global__ void  LPlusTimes_ZGD(double *rhs, double * __restrict__ phi_out, double * __restrict__ ell_plus,
                                int num_zones, int num_groups, int num_local_directions, int num_local_groups, int nidx);



__global__ void  sweep_over_hyperplane_ZGD(int sliceID, int *offset, int *ii_jj_kk_z_idx, int num_directions, int num_groups,
                    int local_imax, int local_jmax, int local_kmax, double *dx, double *dy, double *dz, double *rhs, double *phi, double *psi, 
                    double *sigt, Directions *direction, 
                    double *i_plane, double *j_plane, double *k_plane);


int cuda_LTimes_ZGD(double *d_phi, double *h_psi, double *d_ell,
                    int num_zones, int num_groups, int num_local_directions, int num_local_groups, int nidx);


int  cuda_LPlusTimes_ZGD(double *d_rhs, double *h_phi_out, double *h_ell_plus,
                    int num_zones, int num_groups, int num_local_directions, int num_local_groups, int nidx);

int cuda_sweep_ZGD( double *rhs, double *phi, double *psi,  double *sigt, Directions *direction,
                    double *i_plane, double *j_plane, double *k_plane,
                    int *ii_jj_kk_z_idx, int *offset, int *d_offset, double *dx, double *dy, double *dz,
                    int num_zones, int num_directions, int num_groups,
                    int local_imax, int local_jmax, int local_kmax, int Nslices);


int cuda_LTimes_ZGD(double *d_phi, double *h_psi, double *d_ell,
                    int num_zones, int num_groups, int num_local_directions, int num_local_groups, int nidx){

  cudaCheckError();

  dim3 threadsPerBlock(32);

  LTimes_ZGD<<<num_zones,threadsPerBlock,nidx*sizeof(double)>>>(d_phi,h_psi,d_ell,num_zones,num_groups,num_local_directions,num_local_groups,nidx);

  cudaCheckError();

  return 0;
}

__global__ void  LTimes_ZGD(double *phi, double * __restrict__ psi, double * __restrict__ ell,
                                int num_zones, int num_groups, int num_local_directions, int num_local_groups, int nidx){


      extern __shared__ double ss_phi[];

      int z = blockIdx.x;
      double *block_phi = &phi[z*num_groups*nidx];
      double *block_psi = &psi[z*num_local_groups*num_local_directions];


      for(int group = 0;group < num_local_groups;++ group){
         

        for(int nm_offset = threadIdx.x; nm_offset < nidx; nm_offset+=blockDim.x)
           ss_phi[nm_offset] = block_phi[nm_offset+nidx*group];
 
        for (int d = 0; d < num_local_directions; d++) {

          double psi_d = block_psi[d+group*num_local_directions];
           
          for(int nm_offset = threadIdx.x; nm_offset < nidx; nm_offset+=blockDim.x)
            ss_phi[nm_offset] += ell[nm_offset + nidx*d] * psi_d;          
        }
        for(int nm_offset = threadIdx.x; nm_offset < nidx; nm_offset+=blockDim.x)
           block_phi[nm_offset+nidx*group] =  ss_phi[nm_offset];

      }

}



int  cuda_LPlusTimes_ZGD(double *d_rhs, double *h_phi_out, double *d_ell_plus,
                    int num_zones, int num_groups, int num_local_directions, int num_local_groups, int nidx){


  cudaEvent_t start,stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  float time_ms, time_s;
  double *d_phi_out;

  #ifdef CU_TIMING
  cudaEventRecord(start);
  #endif

  cudaMalloc((void **) &d_phi_out,num_zones*nidx * num_groups * sizeof(double));
//  cudaMalloc((void **) &d_ell_plus, nidx * num_local_directions * sizeof(double));
  cudaMemcpy(d_phi_out, h_phi_out,num_zones * nidx * num_groups * sizeof(double), cudaMemcpyHostToDevice);
//  cudaCheckError();
//  cudaMemcpy(d_ell_plus, h_ell_plus, nidx * num_local_directions * sizeof(double), cudaMemcpyHostToDevice);
  cudaCheckError();

  #ifdef CU_TIMING
  cudaEventRecord(stop);
  cudaDeviceSynchronize();
  cudaCheckError();
  cudaEventElapsedTime(&time_ms,start,stop);
  time_s=time_ms*.001;
  printf("ZGD: time to copy d_phi_out+d_ell_plus H2D: %g [s]\n",time_s);
  #endif
  #ifdef CU_TIMING
  cudaEventRecord(start);
  #endif

  dim3 threadsPerBlock(32);

  LPlusTimes_ZGD<<<num_zones,threadsPerBlock,num_local_directions*sizeof(double)>>>(d_rhs,d_phi_out,d_ell_plus,num_zones,num_groups,num_local_directions,num_local_groups,nidx);

  #ifdef CU_TIMING
  cudaEventRecord(stop);
  cudaDeviceSynchronize();
  cudaCheckError();
  cudaEventElapsedTime(&time_ms,start,stop);
  time_s=time_ms*.001;
  printf("ZGD: time to LPlusTimes_ZGD (GPU): %g [s]\n",time_s);
  #endif


  cudaFree(d_phi_out);
//  cudaFree(d_ell_plus);

  return 0;

}



__global__ void  LPlusTimes_ZGD(double *rhs, double * __restrict__ phi_out,
                                double * __restrict__ ell_plus,
                                int num_zones, int num_groups, int num_local_directions, int num_local_groups, int nidx){


      int z = blockIdx.x;
      extern __shared__ double rhs_acc[];
      double *block_rhs =  &rhs[z*num_groups*num_local_directions];
      double *block_phi_out = &phi_out[z*num_groups*nidx];

      for(int group = 0; group < num_local_groups;++ group){
        for (int d = threadIdx.x; d < num_local_directions; d+=blockDim.x) {
          rhs_acc[d] = 0.0;

          for(int nm_offset = 0;nm_offset < nidx;++nm_offset)
             rhs_acc[d] += ell_plus[nm_offset + d*nidx] * block_phi_out[nm_offset + group*nidx];
          
          block_rhs[d+num_local_directions*group] += rhs_acc[d];
        }
      }
}




int cuda_sweep_ZGD( double *d_rhs, double *h_phi, double *h_psi, double *d_sigt, Directions *d_direction,
                    double *h_i_plane, double *h_j_plane, double *h_k_plane,
                    int *d_ii_jj_kk_z_idx, int *h_offset,  int *d_offset, double *d_dx, double *d_dy, double *d_dz,
                    int num_zones, int num_directions, int num_groups,
                    int local_imax, int local_jmax, int local_kmax, int Nslices){


  size_t N;
  size_t groups_dirs;
  double *d_phi;
  float time_ms, time_s;
  

  N = num_zones * num_directions * num_groups;
  groups_dirs = num_directions * num_groups;


  cudaEvent_t start,stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);



#ifdef USE_PSI_HOST_MEM
  double *d_psi = h_psi;
#else
  double *d_psi;
  #ifdef CU_TIMING
  cudaEventRecord(start);
  #endif
  cudaMalloc((void **) &d_psi, N*sizeof(double));
  cudaMemcpy(d_psi, h_psi,   N*sizeof(double), cudaMemcpyHostToDevice);
  cudaCheckError();

  #ifdef CU_TIMING
  cudaEventRecord(stop);
  cudaDeviceSynchronize();
  cudaCheckError();
  cudaEventElapsedTime(&time_ms,start,stop);
  time_s=time_ms*.001;
  printf("ZGD: time to copy PSI H2D: %g [s]\n",time_s);
  #endif
#endif



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




//  cudaMalloc((void **) &d_phi, N*sizeof(double));

//  cudaMemcpy(d_phi, h_phi,   N*sizeof(double), cudaMemcpyHostToDevice);


  cudaCheckError();
 

  cudaFuncSetCacheConfig(sweep_over_hyperplane_ZGD, cudaFuncCachePreferL1);

//call cuda kernel to sweep over hyperplanes(slices)

  dim3 threadsPerBlock(32,4);

  for (int slice = 0; slice < Nslices; slice++){
    
     #ifdef CU_TIMING
     cudaEventRecord(start);                                             
     #endif

     dim3 numBlocks = h_offset[slice+1] - h_offset[slice];
     sweep_over_hyperplane_ZGD<<<numBlocks,threadsPerBlock, num_groups*sizeof(double) >>>(slice,d_offset,d_ii_jj_kk_z_idx,num_directions,num_groups,local_imax,local_jmax,local_kmax,
                                                          d_dx, d_dy, d_dz, d_rhs, d_phi, d_psi, d_sigt, d_direction, d_i_plane, d_j_plane, d_k_plane);
     #ifdef CU_TIMING
     cudaEventRecord(stop);                                              
     cudaDeviceSynchronize();                                            
     cudaCheckError();                                                   
     float time_ms, time_s;                                              
     cudaEventElapsedTime(&time_ms,start,stop);                          
     time_s=time_ms*.001;                                                
     printf("ZGD: #blocks=%d, time=%g [s]\n",h_offset[slice+1] - h_offset[slice],time_s);
     #endif
     
  }

  #ifdef CU_TIMING
  cudaEventRecord(start);
  #endif

#ifndef USE_IJK_PLANE_HOST_MEM
    cudaMemcpy(h_i_plane, d_i_plane, i_plane_zones * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_j_plane, d_j_plane, j_plane_zones * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_k_plane, d_k_plane, k_plane_zones * sizeof(double), cudaMemcpyDeviceToHost);
#endif

  #ifdef CU_TIMING
  cudaEventRecord(stop);
  cudaDeviceSynchronize();
  cudaCheckError();
  cudaEventElapsedTime(&time_ms,start,stop);
  time_s=time_ms*.001;
  printf("ZGD: time to copy ijk_PLANE D2H: %g [s]\n",time_s);
  #endif


#ifndef USE_PSI_HOST_MEM

#ifdef CU_TIMING
  cudaEventRecord(start);
#endif

  cudaMemcpy(h_psi,     d_psi, N*sizeof(double),                   cudaMemcpyDeviceToHost);

#ifdef CU_TIMING
  cudaEventRecord(stop);
  cudaDeviceSynchronize();
  cudaCheckError();
  cudaEventElapsedTime(&time_ms,start,stop);
  time_s=time_ms*.001;
  printf("ZGD: time to copy PSI D2H: %g [s]\n",time_s);
#endif

  cudaFree(d_psi);
#endif


//  cudaFree(d_phi);

#ifndef USE_IJK_PLANE_HOST_MEM
  cudaFree(d_i_plane);
  cudaFree(d_j_plane);
  cudaFree(d_k_plane);
#endif

  cudaCheckError();

  return 0;
}

#define I_PLANE_INDEX(j, k) ((k)*(local_jmax) + (j))
#define J_PLANE_INDEX(i, k) ((k)*(local_imax) + (i))
#define K_PLANE_INDEX(i, j) ((j)*(local_imax) + (i))


__global__ void sweep_over_hyperplane_ZGD(int sliceID, int *offset, int *ii_jj_kk_z_idx, int num_directions, int num_groups,
                    int local_imax, int local_jmax, int local_kmax, double *dx, double *dy, double *dz, double *rhs, double *phi, double *psi, 
                    double *sigt, Directions *direction,
                    double *i_plane, double *j_plane, double *k_plane){


 
//each block will process one element 
      int element = offset[sliceID] + blockIdx.x;
      if (element > offset[sliceID+1]) return; 

      int i    = ii_jj_kk_z_idx[element*4];
      int j    = ii_jj_kk_z_idx[element*4+1];
      int k    = ii_jj_kk_z_idx[element*4+2];
      int z    = ii_jj_kk_z_idx[element*4+3];
      int dir_grp = num_directions*num_groups;
      int I_P_I = I_PLANE_INDEX(j, k);
      int J_P_I = J_PLANE_INDEX(i, k);
      int K_P_I = K_PLANE_INDEX(i, j);
      double two_inv_dxi = 2.0/dx[i + 1];
      double two_inv_dyj = 2.0/dy[j + 1];
      double two_inv_dzk = 2.0/dz[k + 1];
     

      double * KRESTRICT  block_rhs = &rhs[z*dir_grp];
//    double * KRESTRICT  block_phi = &phi[z*num_directions*num_groups];
      double * KRESTRICT  block_psi = &psi[z*dir_grp];
      double * KRESTRICT  block_sigt = &sigt[z*num_groups];

      double * KRESTRICT psi_lf_z = &i_plane[I_P_I*dir_grp]; 
      double * KRESTRICT psi_fr_z = &j_plane[J_P_I*dir_grp]; 
      double * KRESTRICT psi_bo_z = &k_plane[K_P_I*dir_grp]; 

//      extern __shared__ double xyzcos_dxyz[];

      
//      for( int d = threadIdx.x + threadIdx.y*blockDim.x;  d < num_directions; d = d + blockDim.x*blockDim.y){//not working - BUG?
//      for( int d = threadIdx.x;  d < num_directions; d = d + blockDim.x){
//         xyzcos_dxyz[3*d]   = direction[d].xcos * two_inv_dxi;
//         xyzcos_dxyz[3*d+1] = direction[d].ycos * two_inv_dyj;
//         xyzcos_dxyz[3*d+2] = direction[d].zcos * two_inv_dzk;
//      }

      for (int group = threadIdx.y; group < num_groups; group += blockDim.y){

          for (int  d = threadIdx.x; d < num_directions; d += blockDim.x){

            int gd = d + group*num_directions;

            double xcos_dxi =  direction[d].xcos * two_inv_dxi; 
            double ycos_dyj =  direction[d].ycos * two_inv_dyj;
            double zcos_dzk =  direction[d].zcos * two_inv_dzk;

            double psi_lf_z_g_d = psi_lf_z[gd];
            double psi_fr_z_g_d = psi_fr_z[gd];
            double psi_bo_z_g_d = psi_bo_z[gd];

            /* Calculate new zonal flux */
            double psi_z_g_d = (block_rhs[gd]
                + psi_lf_z_g_d * xcos_dxi
                + psi_fr_z_g_d * ycos_dyj
                + psi_bo_z_g_d * zcos_dzk)
                / (xcos_dxi + ycos_dyj + zcos_dzk + block_sigt[group]);

            block_psi[gd] = psi_z_g_d;


            /* Apply diamond-difference relationships */
            psi_lf_z[gd] = 2.0 * psi_z_g_d - psi_lf_z_g_d;
            psi_fr_z[gd] = 2.0 * psi_z_g_d - psi_fr_z_g_d;
            psi_bo_z[gd] = 2.0 * psi_z_g_d - psi_bo_z_g_d;
          }
        }
}

