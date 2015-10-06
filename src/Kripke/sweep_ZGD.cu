
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



__global__ void  LTimes_ZGD(double *phi, double * __restrict__ psi, const double * __restrict__ ell,
                                int num_zones, int num_groups, int num_local_directions, int num_local_groups, 
                                int nidx, int group0);

__global__ void  LPlusTimes_ZGD(double *rhs, double * __restrict__ phi_out, double * __restrict__ ell_plus,
                                int num_zones, int num_groups, int num_local_directions, int num_local_groups,
                                int nidx, int group0);

__global__ void scattering_ZGD(int * __restrict__ d_mixed_to_zones, int * __restrict__ d_mixed_material,
                               double * __restrict__ d_mixed_fraction, int * __restrict__ d_mixed_offset,
                               double * __restrict__ d_phi, double *d_phi_out, double * __restrict__ d_sigs0,
                              double * __restrict__ d_sigs1, double * __restrict__ d_sigs2,
                               int * __restrict__ moment_to_coeff, int num_mixed, int num_moments, int num_groups, int num_coeff);

__global__ void scattering_ZGD_step2(const int * __restrict__ mixed_to_zones, const int * __restrict__ mixed_material,
                               const double * __restrict__ mixed_fraction, const int * __restrict__ mixed_offset,
                               const double * __restrict__ phi, double *phi_out, double * __restrict__ sigs0,
                               const double * __restrict__ sigs1, const double * __restrict__ sigs2,
                               const int * __restrict__ moment_to_coeff, int num_mixed, int num_moments, int num_groups, int num_coeff);


__global__ void scattering_ZGD_step3(const int * __restrict__ mixed_to_zones, const int * __restrict__ mixed_material,
                               const double * __restrict__ mixed_fraction, const int * __restrict__ mixed_offset,
                               const double * __restrict__ phi, double *phi_out, double * __restrict__ sigs0,
                               const double * __restrict__ sigs1, const double * __restrict__ sigs2,
                               const int * __restrict__ moment_to_coeff, int num_mixed, int num_moments, int num_groups, int num_coeff);


__global__ void source_ZGD(int * __restrict__ mixed_to_zones, int * __restrict__ mixed_material,
                               double * __restrict__ mixed_fraction, int * __restrict__ mixed_offset,
                               double *phi_out, int num_moments, int num_groups);

__global__ void  sweep_over_hyperplane_ZGD(int sliceID, int * __restrict__ offset, int * __restrict__ ii_jj_kk_z_idx, int num_directions, int num_groups,
                    int local_imax, int local_jmax, int local_kmax, 
                    double * __restrict__ dx, double * __restrict__ dy, double * __restrict__ dz, 
                    double * __restrict__ rhs, double * __restrict__ phi, double * __restrict__ psi, 
                    const double * __restrict__ sigt, const Directions * __restrict__ direction, 
                    double *i_plane, double *j_plane, double *k_plane);


__global__ void LPlusTimes_sweep_over_hyperplane_ZGD(int sliceID, int * __restrict__ offset, int * __restrict__ ii_jj_kk_z_idx, int num_directions, int num_groups, int num_local_groups, int nidx, int group0,
                    int local_imax, int local_jmax, int local_kmax, double * __restrict__ dx, double * __restrict__ dy,
                    double * __restrict__ dz, double *__restrict__ phi_out, double * __restrict__ ell_plus, double * __restrict__ psi,
                    const double * __restrict__ sigt, const Directions * __restrict__ direction,
                    double *i_plane, double *j_plane, double *k_plane);
					
					
					
int cuda_LTimes_ZGD(double *d_phi, double *h_psi, double *d_ell,
                    int num_zones, int num_groups, int num_local_directions, int num_local_groups, 
                    int nidx, int group0);


int  cuda_LPlusTimes_ZGD(double *d_rhs, double *h_phi_out, double *h_ell_plus,
                    int num_zones, int num_groups, int num_local_directions, int num_local_groups, 
                    int nidx, int group0);

int  cuda_scattering_ZGD(int *d_mixed_to_zones, int *d_mixed_material, double *d_mixed_fraction, int *d_mixed_offset,
                      double *d_phi, double *d_phi_out, double *d_sigs0,
                      double *d_sigs1, double *d_sigs2,
                      int *moment_to_coeff, int num_mixed, int num_moments, int num_groups, int num_coeff);

int  cuda_source_ZGD(int *d_mixed_to_zones, int *d_mixed_material, double *d_mixed_fraction, int *d_mixed_offset, 
                     double *d_phi_out, int num_moments,  int num_groups);


int cuda_sweep_ZGD( double *rhs, double *phi, double *psi,  double *sigt, Directions *direction,
                    double *i_plane, double *j_plane, double *k_plane,
                    int *ii_jj_kk_z_idx, int *offset, int *d_offset, double *dx, double *dy, double *dz,
                    int num_zones, int num_directions, int num_groups,
                    int local_imax, int local_jmax, int local_kmax, int Nslices);
int cuda_LPlusTimes_sweep_ZGD( double *phi_out, double *ell_plus,
                    double *psi, double *sigt,  Directions *direction,
                    double *i_plane, double *j_plane, double *k_plane,
                    int *ii_jj_kk_z_idx, int *h_offset, int *d_offset, double *dx, double *dy, double *dz,
                    int num_zones, int num_directions, int num_groups, int num_local_groups, int nidx,
                    int local_imax, int local_jmax, int local_kmax, int Nslices);


int cuda_LTimes_ZGD(double *d_phi, double *h_psi, double *d_ell,
                    int num_zones, int num_groups, int num_local_directions, int num_local_groups, int nidx, int group0){

  cudaCheckError();
  int dim_y = 4;  
  dim3 threadsPerBlock(32,dim_y);

  LTimes_ZGD<<<num_zones,threadsPerBlock,nidx*dim_y*sizeof(double)>>>(d_phi,h_psi,d_ell,num_zones,num_groups,num_local_directions,num_local_groups,nidx,group0);

  cudaCheckError();

  return 0;
}

__global__ void  LTimes_ZGD(double *phi, double * __restrict__ psi, const double * __restrict__ ell,
                                int num_zones, int num_groups, int num_local_directions, int num_local_groups, 
                                int nidx, int group0){


      extern __shared__ double ss_phi[];

      int z = blockIdx.x;
      double * __restrict__ block_phi = &phi[z*num_groups*nidx + group0*nidx];
      double * __restrict__ block_psi = &psi[z*num_local_groups*num_local_directions];

#if 1
     for(int group = threadIdx.y ;group < num_local_groups; group += blockDim.y){

        double *ss_phi_group = &ss_phi[nidx*threadIdx.y];


//   offset = z * zones*groups +
//            g * directions +
//            d;

        for(int nm_offset = threadIdx.x; nm_offset < nidx; nm_offset+=blockDim.x)
           ss_phi_group[nm_offset] = block_phi[nm_offset+nidx*group];

        for (int d = 0; d < num_local_directions; d++) {

          double psi_d = block_psi[d+group*num_local_directions];

          for(int nm_offset = threadIdx.x; nm_offset < nidx; nm_offset+=blockDim.x)
            ss_phi_group[nm_offset] += ell[nm_offset + nidx*d] * psi_d;
        }

        for(int nm_offset = threadIdx.x; nm_offset < nidx; nm_offset+=blockDim.x)
           block_phi[nm_offset+nidx*group] =  ss_phi_group[nm_offset];
      }

#else

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
#endif

}


/*******************/
int  cuda_scattering_ZGD(int *d_mixed_to_zones, int *d_mixed_material, double *d_mixed_fraction, int *d_mixed_offset,
                      double *d_phi, double *d_phi_out, double *d_sigs0, double *d_sigs1, double *d_sigs2,
                      int *d_moment_to_coeff, int num_mixed, int num_moments, int num_groups, int num_coeff){


     int y_dim = 6;
     dim3 threadsPerBlock(32,y_dim);


//     scattering_ZGD<<<480,threadsPerBlock>>>(d_mixed_to_zones,d_mixed_material,d_mixed_fraction,d_mixed_offset,
//                                d_phi,d_phi_out,d_sigs0, d_sigs1, d_sigs2, d_moment_to_coeff,num_mixed,num_moments,num_groups,num_coeff);


     scattering_ZGD_step3<<<480,threadsPerBlock,num_groups*y_dim*sizeof(double)>>>(d_mixed_to_zones,d_mixed_material,d_mixed_fraction,d_mixed_offset,
                                d_phi,d_phi_out,d_sigs0, d_sigs1, d_sigs2, d_moment_to_coeff,num_mixed,num_moments,num_groups,num_coeff);

     cudaCheckError();

    return 0;
}

/*******************/



__global__ void scattering_ZGD(int * __restrict__ mixed_to_zones, int * __restrict__ mixed_material,
                               double * __restrict__ mixed_fraction, int * __restrict__ mixed_offset,
                               double * __restrict__ phi, double *phi_out, double * __restrict__ sigs0,
                               double * __restrict__ sigs1, double * __restrict__ sigs2,
                               int * __restrict__ moment_to_coeff, int num_mixed, int num_moments, int num_groups, int num_coeff){


   int mix_min = mixed_offset[blockIdx.x];
   int mix_max = mixed_offset[blockIdx.x+1];
   double *d_sigs[3];
   d_sigs[0] = sigs0;
   d_sigs[1] = sigs1;
   d_sigs[2] = sigs2;

   for(int mix = mix_min;mix < mix_max; ++mix){
      int zone = mixed_to_zones[mix];
      int material = mixed_material[mix];
      double fraction = mixed_fraction[mix];
      double *sigs_g_gp = d_sigs[material];
      double *phi_z_g = &phi[zone*num_groups*num_moments];
      double *phi_out_z_gp = &phi_out[zone*num_groups*num_moments];

//LG unroll the outer loop
      for(int g = 0; g < num_groups;++g){

        for(int gp = threadIdx.y; gp < num_groups; gp += blockDim.y){

          for(int nm = threadIdx.x; nm < num_moments; nm += blockDim.x){
            // map nm to n
            int n = moment_to_coeff[nm];

            phi_out_z_gp[nm+gp*num_moments] += 
                sigs_g_gp[n + gp*num_coeff + g*num_groups*num_coeff] * 
                phi_z_g[nm + g*num_moments] * fraction;
          }
        }
      }

   }
}



__global__ void scattering_ZGD_step2(const int * __restrict__ mixed_to_zones, const int * __restrict__ mixed_material,
                               const double * __restrict__ mixed_fraction, const int * __restrict__ mixed_offset,
                               const double * __restrict__ phi, double *phi_out, double * __restrict__ sigs0,
                               const double * __restrict__ sigs1, const double * __restrict__ sigs2,
                               const int * __restrict__ moment_to_coeff, int num_mixed, int num_moments, int num_groups, int num_coeff){

   extern __shared__ double phi_z_g_ss[];

   int mix_min = mixed_offset[blockIdx.x];
   int mix_max = mixed_offset[blockIdx.x+1];
   int gid = threadIdx.x + threadIdx.y*blockDim.x;
   const double *d_sigs[3];
   d_sigs[0] = sigs0;
   d_sigs[1] = sigs1;
   d_sigs[2] = sigs2;

   for(int mix = mix_min;mix < mix_max; ++mix){
      int zone = mixed_to_zones[mix];
      int material = mixed_material[mix];
      double fraction = mixed_fraction[mix];
      const double * __restrict__ sigs_g_gp = d_sigs[material];
      const double * __restrict__ phi_z_g = &phi[zone*num_groups*num_moments];
      double * __restrict__ phi_out_z_gp = &phi_out[zone*num_groups*num_moments];

//LG unroll the outer loop
      for(int g = 0; g < num_groups;++g){

        __syncthreads();
        for(int nm = gid; nm < num_moments; nm += blockDim.x*blockDim.y)
          phi_z_g_ss[nm] =  phi_z_g[nm + g*num_moments]*fraction;
         __syncthreads(); 

        for(int gp = threadIdx.y; gp < num_groups; gp += blockDim.y){

          for(int nm = threadIdx.x; nm < num_moments; nm += blockDim.x){
            // map nm to n
            const int n = moment_to_coeff[nm];

            phi_out_z_gp[nm+gp*num_moments] +=
                sigs_g_gp[n + gp*num_coeff + g*num_groups*num_coeff] *
                phi_z_g_ss[nm];
          }
        }
      }

   }
}


__global__ void scattering_ZGD_step3(const int * __restrict__ mixed_to_zones, const int * __restrict__ mixed_material,
                               const double * __restrict__ mixed_fraction, const int * __restrict__ mixed_offset,
                               const double * __restrict__ phi, double *phi_out, double * __restrict__ sigs0,
                               const double * __restrict__ sigs1, const double * __restrict__ sigs2,
                               const int * __restrict__ moment_to_coeff, int num_mixed, int num_moments, int num_groups, int num_coeff){

   extern __shared__ double phi_out_ss[];

   int mix_min = mixed_offset[blockIdx.x];
   int mix_max = mixed_offset[blockIdx.x+1];
   const double *d_sigs[3];
   d_sigs[0] = sigs0;
   d_sigs[1] = sigs1;
   d_sigs[2] = sigs2;

   for(int mix = mix_min;mix < mix_max; ++mix){
      int zone = mixed_to_zones[mix];
      int material = mixed_material[mix];
      double fraction = mixed_fraction[mix];
      const double * __restrict__ sigs_g_gp = d_sigs[material];
      const double * __restrict__ phi_z_g = &phi[zone*num_groups*num_moments];
      double * __restrict__ phi_out_z_gp = &phi_out[zone*num_groups*num_moments];

//LG unroll the outer loop

      for(int gp = threadIdx.y; gp < num_groups; gp += blockDim.y){

        double *phi_out_ss_gp = &phi_out_ss[num_groups*threadIdx.y];
        for(int nm = threadIdx.x; nm < num_moments; nm += blockDim.x)
           phi_out_ss_gp[nm] =  phi_out_z_gp[nm+gp*num_moments];

        for(int g = 0; g < num_groups;++g){

          const int nm_shift = g*num_groups*num_coeff + gp*num_coeff ;
          for(int nm = threadIdx.x; nm < num_moments; nm += blockDim.x){
            // map nm to n
            const int n = moment_to_coeff[nm];

            phi_out_ss_gp[nm] +=
                sigs_g_gp[n + nm_shift] *
                phi_z_g[nm + g*num_moments]*fraction;

          }
        }
        for(int nm = threadIdx.x; nm < num_moments; nm += blockDim.x)
           phi_out_z_gp[nm+gp*num_moments] = phi_out_ss_gp[nm];

      }

   }
}







int  cuda_source_ZGD(int *d_mixed_to_zones, int *d_mixed_material, double *d_mixed_fraction, int *d_mixed_offset,
                      double *d_phi_out, int num_moments, int num_groups){

     dim3 threadsPerBlock(32,1);
     source_ZGD<<<480,threadsPerBlock>>>(d_mixed_to_zones,d_mixed_material,d_mixed_fraction,d_mixed_offset,
                                d_phi_out,num_moments,num_groups);

     cudaCheckError();
    return 0;
}

__global__ void source_ZGD(int * __restrict__ mixed_to_zones, int * __restrict__ mixed_material,
                               double * __restrict__ mixed_fraction, int * __restrict__ mixed_offset,
                               double *phi_out, int num_moments, int num_groups){


   int mix_min = mixed_offset[blockIdx.x];
   int mix_max = mixed_offset[blockIdx.x+1];


   for(int mix = mix_min;mix < mix_max; ++mix){
  
      int material = mixed_material[mix];

      if(material == 0){
          int zone = mixed_to_zones[mix];
          double fraction = mixed_fraction[mix];
          double *phi_out_z = &phi_out[zone*num_moments*num_groups];
          for(int g = threadIdx.x; g < num_groups; g += blockDim.x){
            phi_out_z[g*num_moments] += 1.0 * fraction;
          }
       }
    }

}


int  cuda_LPlusTimes_ZGD(double *d_rhs, double *d_phi_out, double *d_ell_plus,
                    int num_zones, int num_groups, int num_local_directions, int num_local_groups,
                    int nidx, int group0){


  #ifdef CU_TIMING
  cudaEvent_t start,stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  float time_ms, time_s;

  cudaEventRecord(start);
  #endif

  dim3 threadsPerBlock(32);

  LPlusTimes_ZGD<<<num_zones,threadsPerBlock,num_local_directions*sizeof(double)>>>(d_rhs,d_phi_out,d_ell_plus,num_zones,num_groups, num_local_directions, num_local_groups, nidx, group0);

  #ifdef CU_TIMING
  cudaEventRecord(stop);
  cudaDeviceSynchronize();
  cudaCheckError();
  cudaEventElapsedTime(&time_ms,start,stop);
  time_s=time_ms*.001;
  printf("ZGD: time to LPlusTimes_ZGD (GPU): %g [s]\n",time_s);
  #endif

  return 0;

}



__global__ void  LPlusTimes_ZGD(double *rhs, double * __restrict__ phi_out,
                                double * __restrict__ ell_plus,
                                int num_zones, int num_groups, int num_local_directions, int num_local_groups, 
                                int nidx, int group0){


//LG consider running  for(int nm_offset  in parallel and using reduction within a warp

      int z = blockIdx.x;
      extern __shared__ double rhs_acc[];
      double *block_rhs =  &rhs[z*num_local_groups*num_local_directions];


//   offset = z * zones*groups +
//            g * directions +
//            d;


      double *block_phi_out = &phi_out[z*num_groups*nidx + group0*nidx];

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

  cudaCheckError();
 

//  cudaFuncSetCacheConfig(sweep_over_hyperplane_ZGD, cudaFuncCachePreferL1);

//call cuda kernel to sweep over hyperplanes(slices)

  dim3 threadsPerBlock(32,8);

  for (int slice = 0; slice < Nslices; slice++){
    
     #ifdef CU_TIMING
     cudaEventRecord(start);                                             
     #endif

     dim3 numBlocks = h_offset[slice+1] - h_offset[slice];
     sweep_over_hyperplane_ZGD<<<numBlocks,threadsPerBlock>>>(slice,d_offset,d_ii_jj_kk_z_idx,num_directions,num_groups,local_imax,local_jmax,local_kmax,
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


__global__ void  sweep_over_hyperplane_ZGD(int sliceID, int * __restrict__ offset, int * __restrict__ ii_jj_kk_z_idx, int num_directions, int num_groups,
                    int local_imax, int local_jmax, int local_kmax,
                    double * __restrict__ dx, double * __restrict__ dy, double * __restrict__ dz,
                    double * __restrict__ rhs, double * __restrict__ phi, double * __restrict__ psi,
                    const double * __restrict__ sigt, const Directions * __restrict__ direction,
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
      const double * KRESTRICT  block_sigt = &sigt[z*num_groups];

      double * KRESTRICT psi_lf_z = &i_plane[I_P_I*dir_grp]; 
      double * KRESTRICT psi_fr_z = &j_plane[J_P_I*dir_grp]; 
      double * KRESTRICT psi_bo_z = &k_plane[K_P_I*dir_grp]; 

     

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

int cuda_LPlusTimes_sweep_ZGD( double *d_phi_out, double *d_ell_plus, double *h_psi, double *d_sigt, Directions *d_direction,
                    double *h_i_plane, double *h_j_plane, double *h_k_plane,
                    int *d_ii_jj_kk_z_idx, int *h_offset,  int *d_offset, double *d_dx, double *d_dy, double *d_dz,
                    int num_zones, int num_directions, int num_groups, int num_local_groups,
                    int nidx, int group0,
                    int local_imax, int local_jmax, int local_kmax, int Nslices){


  size_t N;
  size_t groups_dirs;
  float time_ms, time_s;
  
  N = num_zones * num_directions * num_local_groups;
  groups_dirs = num_directions * num_local_groups;

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

  cudaCheckError();
 

  cudaFuncSetCacheConfig(sweep_over_hyperplane_ZGD, cudaFuncCachePreferL1);

//call cuda kernel to sweep over hyperplanes(slices)
  int dim_y = 8; 
  dim3 threadsPerBlock(32,dim_y);

  for (int slice = 0; slice < Nslices; slice++){
    
     #ifdef CU_TIMING
     cudaEventRecord(start);                                             
     #endif

     dim3 numBlocks = h_offset[slice+1] - h_offset[slice];
     LPlusTimes_sweep_over_hyperplane_ZGD<<<numBlocks,threadsPerBlock,num_directions*dim_y*sizeof(double)>>>(slice,d_offset,d_ii_jj_kk_z_idx,num_directions,num_groups,num_local_groups,nidx,group0, local_imax,local_jmax,local_kmax,
                                                          d_dx, d_dy, d_dz, d_phi_out, d_ell_plus, d_psi, d_sigt, d_direction, d_i_plane, d_j_plane, d_k_plane);
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


__global__ void LPlusTimes_sweep_over_hyperplane_ZGD(int sliceID, int * __restrict__ offset, int * __restrict__ ii_jj_kk_z_idx, int num_directions, int num_groups, int num_local_groups, int nidx, int group0,
                    int local_imax, int local_jmax, int local_kmax, double * __restrict__ dx, double * __restrict__ dy, 
                    double * __restrict__ dz, double *__restrict__ phi_out, double * __restrict__ ell_plus, double * __restrict__ psi, 
                    const double * __restrict__ sigt, const Directions * __restrict__ direction,
                    double *i_plane, double *j_plane, double *k_plane){
 
     extern __shared__ double rhs_local[];
 
//each block will process one element 
      int element = offset[sliceID] + blockIdx.x;
      if (element > offset[sliceID+1]) return; 

      int i    = ii_jj_kk_z_idx[element*4];
      int j    = ii_jj_kk_z_idx[element*4+1];
      int k    = ii_jj_kk_z_idx[element*4+2];
      int z    = ii_jj_kk_z_idx[element*4+3];
      int dir_grp = num_directions*num_local_groups;
      int I_P_I = I_PLANE_INDEX(j, k);
      int J_P_I = J_PLANE_INDEX(i, k);
      int K_P_I = K_PLANE_INDEX(i, j);
      double two_inv_dxi = 2.0/dx[i + 1];
      double two_inv_dyj = 2.0/dy[j + 1];
      double two_inv_dzk = 2.0/dz[k + 1];
     
      double * KRESTRICT  block_psi = &psi[z*dir_grp];
      const double * KRESTRICT  block_sigt = &sigt[z*num_local_groups];
      double * KRESTRICT psi_lf_z = &i_plane[I_P_I*dir_grp]; 
      double * KRESTRICT psi_fr_z = &j_plane[J_P_I*dir_grp]; 
      double * KRESTRICT psi_bo_z = &k_plane[K_P_I*dir_grp];
 
      const double * KRESTRICT block_phi_out = &phi_out[z*num_groups*nidx + group0*nidx];

      for (int group = threadIdx.y; group < num_local_groups; group += blockDim.y){

	  double * __restrict__ rhs_local_group = &rhs_local[threadIdx.y * num_directions];
          for (int d = threadIdx.x; d < num_directions; d+=blockDim.x) {
            double sum  = 0.0;
            for(int nm_offset = 0;nm_offset < nidx;++nm_offset)
               sum += ell_plus[nm_offset + d*nidx] * block_phi_out[nm_offset + group*nidx];
            rhs_local_group[d] = sum;     
          }	  
	 
          for (int  d = threadIdx.x; d < num_directions; d += blockDim.x){

            int gd = d + group*num_directions;

            double xcos_dxi =  direction[d].xcos * two_inv_dxi; 
            double ycos_dyj =  direction[d].ycos * two_inv_dyj;
            double zcos_dzk =  direction[d].zcos * two_inv_dzk;

            double psi_lf_z_g_d = psi_lf_z[gd];
            double psi_fr_z_g_d = psi_fr_z[gd];
            double psi_bo_z_g_d = psi_bo_z[gd];

            /* Calculate new zonal flux */
            double psi_z_g_d = (rhs_local_group[d]
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
