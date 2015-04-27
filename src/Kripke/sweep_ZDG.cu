
#include <stdio.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "Directions.h"


#define KRESTRICT __restrict__

#define USE_PSI_HOST_MEM


//#define USE_IJK_PLANE_HOST_MEM



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

__global__ void  LTimes_ZDG(double *phi, double * __restrict__ psi, double * __restrict__ ell,
                                int num_zones, int num_groups, int num_local_directions, int num_local_groups, int nidx);

__global__ void  LPlusTimes_ZDG(double *rhs, double * __restrict__ phi_out, double * __restrict__ ell_plus,
                                int num_zones, int num_groups, int num_local_directions, int num_local_groups, int nidx);


__global__ void  sweep_over_hyperplane_ZDG(int const sliceID, int * __restrict__ offset,  int * __restrict__ ii_jj_kk_z_idx, int const num_directions, int const num_groups,
                    int const local_imax, int const local_jmax, int const local_kmax, double * __restrict__ dx, double * __restrict__ dy, double * __restrict__ dz, 
                    double * __restrict__ rhs, double *phi, double *psi, 
                    double * __restrict__ sigt, Directions * __restrict__ direction, 
                    double *i_plane, double *j_plane, double *k_plane);


int cuda_LTimes_ZDG(double *d_phi, double *h_psi, double *d_ell,
                    int num_zones, int num_groups, int num_local_directions, int num_local_groups, int nidx);

int  cuda_LPlusTimes_ZDG(double *rhs, double *phi_out, double *ell_plus,
                    int num_zones, int num_groups, int num_local_directions, int num_local_groups, int nidx);

int cuda_sweep_ZDG( double *rhs, double *phi, double *psi,  double *sigt, Directions *direction,
                    double *i_plane, double *j_plane, double *k_plane,
                    int *ii_jj_kk_z_idx, int *offset, int *d_offset, double *dx, double *dy, double *dz,
                    int num_zones, int num_directions, int num_groups,
                    int local_imax, int local_jmax, int local_kmax, int Nslices);


/*******************/

int cuda_LTimes_ZDG(double *d_phi, double *h_psi, double *d_ell,
                    int num_zones, int num_groups, int num_local_directions, int num_local_groups, int nidx){

  
  cudaCheckError();

  dim3 threadsPerBlock(32);

  LTimes_ZDG<<<num_zones,threadsPerBlock,(num_local_groups+nidx)*sizeof(double)>>>(d_phi,h_psi,d_ell,num_zones,num_groups,num_local_directions,num_local_groups,nidx);

  cudaCheckError();


  return 0;
}


__global__ void  LTimes_ZDG(double *phi, double * __restrict__ psi, double * __restrict__ ell,
                                int num_zones, int num_groups, int num_local_directions, int num_local_groups, int nidx){

      extern __shared__  double tmp_psi[];

//need more shared memory
      double *tmp_ell_d_nm = &tmp_psi[num_local_groups];


      int z = blockIdx.x;
      double *block_phi = &phi[z*num_groups*nidx];
      double *block_psi = &psi[z*num_local_groups*num_local_directions];


      for (int d = 0; d < num_local_directions; d++) {
        
        for (int group = threadIdx.x ; group < num_local_groups; group+=blockDim.x )
           tmp_psi[group] = block_psi[group + d*num_local_groups];

//use shared memory for tmp_ell_d_nm
        for (int nm_offset = threadIdx.x ; nm_offset < nidx; nm_offset += blockDim.x )
           tmp_ell_d_nm[nm_offset] =  ell[nm_offset+d*nidx];
        __syncthreads();

        for(int nm_offset = 0;nm_offset < nidx;++nm_offset){
          double ell_d_nm =  tmp_ell_d_nm[nm_offset];//  ell[nm_offset+d*nidx]; 

          for (int group = threadIdx.x ; group < num_local_groups; group+=blockDim.x )
            block_phi[group + num_groups*nm_offset] += ell_d_nm * tmp_psi[group];
                   
        }
        
    }
}

/*******************/

int  cuda_LPlusTimes_ZDG(double *d_rhs, double *h_phi_out, double *d_ell_plus,       
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
  cudaMemcpy(d_phi_out, h_phi_out,num_zones * nidx * num_groups * sizeof(double), cudaMemcpyHostToDevice);
  cudaCheckError();

  #ifdef CU_TIMING
  cudaEventRecord(stop);
  cudaDeviceSynchronize();
  cudaCheckError();
  cudaEventElapsedTime(&time_ms,start,stop);
  time_s=time_ms*.001;
  printf("ZDG: time to copy d_phi_out+d_ell_plus H2D: %g [s]\n",time_s);
  #endif


  #ifdef CU_TIMING
  cudaEventRecord(start);
  #endif

  dim3 threadsPerBlock(32);

  LPlusTimes_ZDG<<<num_zones,threadsPerBlock,num_local_groups*sizeof(double)>>>(d_rhs,d_phi_out,d_ell_plus,num_zones,num_groups,num_local_directions,num_local_groups,nidx);

  #ifdef CU_TIMING
  cudaEventRecord(stop);
  cudaDeviceSynchronize();
  cudaCheckError();
  cudaEventElapsedTime(&time_ms,start,stop);
  time_s=time_ms*.001;
  printf("ZDG: time to LPlusTimes_ZDG (GPU): %g [s]\n",time_s);
  #endif


  cudaFree(d_phi_out);
  //cudaFree(d_ell_plus);

  return 0;

}



__global__ void  LPlusTimes_ZDG(double *rhs, double * __restrict__ phi_out,
                                double * __restrict__ ell_plus,
                                int num_zones, int num_groups, int num_local_directions, int num_local_groups, int nidx){

      extern __shared__  double tmp_rhs[];

      int z = blockIdx.x;
      double *block_rhs =  &rhs[z*num_groups*num_local_directions];
      double *block_phi_out = &phi_out[z*num_groups*nidx];

      for (int d = 0; d < num_local_directions; d++) {

        for (int group = threadIdx.x ; group < num_local_groups; group+=blockDim.x )
           tmp_rhs[group] = 0.0;

        //unroll if possible 
        for(int nm_offset = 0;nm_offset < nidx;++nm_offset){
          double ell_plus_d_n_m = ell_plus[nm_offset + d * nidx];
          for (int group = threadIdx.x ; group < num_local_groups; group+=blockDim.x ) 
            tmp_rhs[group] += ell_plus_d_n_m * block_phi_out[group + num_groups*nm_offset];
        }
        //not sure if __syncthreads  is needed, probably not;  
        __syncthreads();
        //copy from shared to global memory 
        for (int group = threadIdx.x ; group < num_local_groups; group+=blockDim.x )
            block_rhs[group + d * num_local_groups] = tmp_rhs[group];
      }
}






int cuda_sweep_ZDG( double *d_rhs, double *h_phi, double *h_psi, double *d_sigt, Directions *d_direction, 
                    double *h_i_plane, double *h_j_plane, double *h_k_plane,
                    int *d_ii_jj_kk_z_idx, int *h_offset,  int *d_offset, double *d_dx, double *d_dy, double *d_dz,
                    int num_zones, int num_directions, int num_groups,
                    int local_imax, int local_jmax, int local_kmax, int Nslices){

  size_t N;
  size_t groups_dirs;
  double *d_phi;
  //double *d_dx, *d_dy, *d_dz;
  float time_ms, time_s;
  static int  INIT_FLAG_IJK_PLANE = 0;


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

  #ifdef CU_TIMING
  cudaEventRecord(stop);
  cudaDeviceSynchronize();
  cudaCheckError();
  cudaEventElapsedTime(&time_ms,start,stop);
  time_s=time_ms*.001;
  printf("ZDG: time to copy PSI H2D: %g [s]\n",time_s);
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
  printf("ZDG: time to copy ijk_PLANE H2D: %g [s]\n",time_s);
  #endif

#endif


//  cudaMemcpy(d_phi, h_phi,   N*sizeof(double), cudaMemcpyHostToDevice);


  cudaFuncSetCacheConfig(sweep_over_hyperplane_ZDG, cudaFuncCachePreferL1);

//call cuda kernel to sweep over hyperplanes(slices)

  dim3 threadsPerBlock(32,12);

  for (int slice = 0; slice < Nslices; slice++){

#ifdef CU_TIMING    
     cudaEventRecord(start); 
#endif                                            
     dim3 numBlocks = h_offset[slice+1] - h_offset[slice];
     sweep_over_hyperplane_ZDG<<<numBlocks,threadsPerBlock>>>(slice,d_offset,d_ii_jj_kk_z_idx,num_directions,num_groups,local_imax,local_jmax,local_kmax,
                                                          d_dx, d_dy, d_dz, d_rhs, d_phi, d_psi, d_sigt, d_direction, d_i_plane, d_j_plane, d_k_plane);

#ifdef CU_TIMING
    cudaEventRecord(stop);                                              
    cudaDeviceSynchronize();                                            
    cudaCheckError();                                                   
    float time_ms, time_s;                                              
    cudaEventElapsedTime(&time_ms,start,stop);                          
    time_s=time_ms*.001;                                                
    printf("ZDG: #blocks=%d, time=%g [s]\n",h_offset[slice+1] - h_offset[slice],time_s);
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
  printf("ZDG: time to copy ijk_PLANE D2H: %g [s]\n",time_s);
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
  printf("ZDG: time to copy PSI D2H: %g [s]\n",time_s);
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


__global__ void sweep_over_hyperplane_ZDG(int const sliceID, int * __restrict__ offset, int * __restrict__ ii_jj_kk_z_idx,  int const num_directions, int const num_groups,
                    int const local_imax, int const local_jmax, int const local_kmax, double * __restrict__ dx, double * __restrict__ dy, double * __restrict__ dz, 
                    double * __restrict__ rhs, double *phi, double *psi, 
                    double * __restrict__ sigt, Directions * __restrict__ direction,
                    double *i_plane, double *j_plane, double *k_plane){

//each block will process one element 
      int element = offset[sliceID] + blockIdx.x;
      if (element > offset[sliceID+1]) return; 
      

      int i    = ii_jj_kk_z_idx[element*4];
      int j    = ii_jj_kk_z_idx[element*4+1];
      int k    = ii_jj_kk_z_idx[element*4+2];
      int z    = ii_jj_kk_z_idx[element*4+3];
      int I_P_I = I_PLANE_INDEX(j, k);
      int J_P_I = J_PLANE_INDEX(i, k);
      int K_P_I = K_PLANE_INDEX(i, j);
      double two_inv_dxi = 2.0/dx[i + 1];
      double two_inv_dyj = 2.0/dy[j + 1];
      double two_inv_dzk = 2.0/dz[k + 1];

      int dir_grp = num_directions*num_groups;

     

      double * KRESTRICT  block_rhs = &rhs[z*dir_grp];
//    double * KRESTRICT  block_phi = &phi[z*num_directions*num_groups];
      double * KRESTRICT  block_psi = &psi[z*dir_grp];
      double * KRESTRICT  block_sigt = &sigt[z*num_groups];

      double * KRESTRICT psi_lf_z_d = &i_plane[I_P_I*dir_grp]; // = i_plane.ptr(0, 0, I_P_I);
      double * KRESTRICT psi_fr_z_d = &j_plane[J_P_I*dir_grp]; // = j_plane.ptr(0, 0, J_P_I);
      double * KRESTRICT psi_bo_z_d = &k_plane[K_P_I*dir_grp]; // = k_plane.ptr(0, 0, K_P_I);


      for (int d = threadIdx.y; d < num_directions; d += blockDim.y){
          for (int group = threadIdx.x; group < num_groups; group += blockDim.x){

            double xcos = direction[d].xcos;
            double ycos = direction[d].ycos;
            double zcos = direction[d].zcos;

            int gd = group + d*num_groups;

            double xcos_dxi = xcos * two_inv_dxi;
            double ycos_dyj = ycos * two_inv_dyj;
            double zcos_dzk = zcos * two_inv_dzk;

            double psi_lf_z_d_group = psi_lf_z_d[gd];
            double psi_fr_z_d_group = psi_fr_z_d[gd];
            double psi_bo_z_d_group = psi_bo_z_d[gd];

            /* Calculate new zonal flux */
            double psi_z_d_g = (block_rhs[gd]
                + psi_lf_z_d_group * xcos_dxi
                + psi_fr_z_d_group * ycos_dyj
                + psi_bo_z_d_group * zcos_dzk)
                / (xcos_dxi + ycos_dyj + zcos_dzk + block_sigt[group]);

            block_psi[gd] = psi_z_d_g;

            /* Apply diamond-difference relationships */
            psi_lf_z_d[gd] = 2.0 * psi_z_d_g - psi_lf_z_d_group;
            psi_fr_z_d[gd] = 2.0 * psi_z_d_g - psi_fr_z_d_group;
            psi_bo_z_d[gd] = 2.0 * psi_z_d_g - psi_bo_z_d_group;
          }
      }
}





