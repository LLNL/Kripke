#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <cuda_runtime.h>
#include <cublas_v2.h>


#define  blk_m 25
#define  blk_n 16 
#define  blk_k 16


/* function declarations */
__global__ void dgemm_custom_simple_1block_batch_kernel( int transa, int transb, 
		  				         int *mlist, int *nlist, int *klist,
						         double alpha, 
						         const double **Alist, int *ldalist,				 
					   	         const double **Blist, int *ldblist,
						         double beta,
			  		                 double **Clist, int *ldclist,
							 int nbatch,
							 int *d_counter 
							 );
			             


extern "C" {
  /* DGEMM wrapper */
void dgemm_custom_simple_1block_batch (cudaStream_t stream, 
  		     	               cublasOperation_t transa, 
				       cublasOperation_t transb, 
				       int *mlist, int *nlist, int *klist, 
				       const double *alpha, 
				       const double **Alist, int *ldalist, 
				       const double **Blist, int *ldblist, 
				       const double *beta, 
				       double **Clist, int *ldclist, int nbatch)                                       
  {

    cudaError_t cuerr = cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeEightByte);
    if ( cuerr ) {
      printf ("error cudaFuncSetSharedMemConfig \n");
      abort();
    }

    /* set transpose */
    int ta = (transa == CUBLAS_OP_T) || (transa == CUBLAS_OP_C);
    int tb = (transb == CUBLAS_OP_T) || (transb == CUBLAS_OP_C);

    /* grid size (batch size) */
    //dim3 grid(nbatch);   
    dim3 grid(120);

    /* create counter on device */
    static int *d_counter = NULL;
    if ( ! d_counter ) {
      cudaMalloc ( (void**)&d_counter, sizeof(int) );
    }
    cudaMemset ( d_counter, 0, sizeof(int) );

    /* block size */
    dim3 threads( 16, 16 );

    /* batched kernel */
    dgemm_custom_simple_1block_batch_kernel<<< grid, threads, 0, stream >>> (ta, tb, mlist, nlist, klist, *alpha, Alist, ldalist, Blist, ldblist, *beta, Clist, ldclist, nbatch, d_counter);
    
    //cudaDeviceSynchronize();
   
  }                                      
}


/* batched DGEMM kernel */
__global__ void dgemm_custom_simple_1block_batch_kernel( int transa, int transb,
		                                         int *mlist, int *nlist, int *klist,
                    			                 double alpha,
		                                         const double **Alist, int *ldalist,
                    			                 const double **Blist, int *ldblist,
		                                         double beta,
                    			                 double **Clist, int *ldclist,
							 int nbatch, 
							 int *d_counter)
{

    /* allocate shared memory */
  __shared__ double As[(blk_m+1)*blk_k];
  __shared__ double Bs[4][(blk_k+1)*blk_n];
  __shared__ int s_count;
  
  /* set thread index */
  int tx = threadIdx.x;
  int ty = threadIdx.y;

  int bdy = blockDim.y;
  int bdx = blockDim.x;

    /* local variables */ 
    int im, in, ik, i, id;
    int colA, colB, colC;
    int rowA, rowB, rowC;
    double ABSUM = 0.0;
    double ABSUM1 = 0.0;
    double ABSUM2 = 0.0;
    double ABSUM3 = 0.0;
 
    if ( tx == 0 && ty == 0 ) {
       s_count = atomicAdd(d_counter,1);
    } 

    __syncthreads();

    int count = s_count;

    if ( s_count >= nbatch ) return;

    /* get dimensions for current block (DGEMM) */
    int m = mlist[count];
    int n = nlist[count]; 
    int k = klist[count];
    int lda = ldalist[count];
    int ldb = ldblist[count];
    int ldc = ldclist[count];
    const double *A = Alist[count];
    const double *B = Blist[count];
    double *C = Clist[count];

    /*load L*/
    for ( ik=0; ik<k; ik+=bdy ) {
      for ( im=0; im<m; im+=bdx ) {
	if ((tx + im < m) && (ty + ik < k)) {
	  As[(tx+im)+(blk_m+1)*(ty+ik)] = __ldg(A+ (tx+im)+lda*(ty+ik) );
	}
      }
    }
    
  int quad = tx/8+2*ty/8;

    while ( s_count < nbatch ) {
    
      /* get dimensions for current block (DGEMM) */
      int m = mlist[count];
      int n = nlist[count]; 
      int k = klist[count];
      int lda = ldalist[count];
      int ldb = ldblist[count];
      int ldc = ldclist[count];
      const double *A = Alist[count];
      const double *B = Blist[count];
      double *C = Clist[count];
      double A0, A1, B0, B1;

      __syncthreads();   
      
      /* early exit */
      if(!m || !n || !k) return;

      /* load B */
      for ( ik=0; ik<n; ik+=bdy ) {
	for ( im=0; im<k; im+=bdx ) {
	  if ((tx + im < k) && (ty + ik < n)) {
	    Bs[quad][(tx+im)+(blk_k+1)*(ty+ik)] = __ldg(B+(tx+im)+ldb*(ty+ik) );
	  }
	}
      }

      /* synchronize threads */
      __syncthreads();
    
      for ( in=0; in<n; in+=2*bdy ) {
	for ( im=0; im<k; im+=2*bdx ) {

	  rowC = tx+im;
	  colC = ty+in;

	  ABSUM = 0.0;
	  ABSUM1 = 0.0;
	  ABSUM2 = 0.0;
	  ABSUM3 = 0.0;

	  for(i = 0; i < k; i++){
	    if ((rowC < m) && (colC < n)) {
	      /* compute A*B */
	      A0= As[rowC+(blk_m+1)*i];
	      B0 = Bs[quad][i+(blk_k+1)*colC];
	      //ABSUM += As[rowC+(blk_m+1)*i]*Bs[0][i+(blk_k+1)*colC];
	      ABSUM += A0*B0;
	    }
	    if ((rowC+bdx < m) && (colC < n)) {
	      /* compute A*B */
	      A1 = As[(rowC+bdx)+(blk_m+1)*i];
	      ABSUM1 += A1*B0;
	    }
	    if ((rowC < m) && (colC+bdy < n)) {
	      /* compute A*B */
	      B1 = Bs[quad][i+(blk_k+1)*(colC+bdy)];
	      ABSUM2 += A0*B1;
	    }
	    if ((rowC+bdx < m) && (colC+bdy < n)) {
	      /* compute A*B */
	      //ABSUM3 += As[(rowC+bdx)+(blk_m+1)*i]*Bs[0][i+(blk_k+1)*(colC+bdy)];
	      ABSUM3 += A1*B1;
	    }
	  }

	  if ((rowC < m) && (colC < n)) {
	    C[ rowC + ldc*colC ] = alpha *ABSUM + beta*C[rowC+ldc*colC];
	  }


	  if ((rowC+bdx < m) && (colC < n)) {
	    C[ (rowC+bdx) + ldc*colC ] = alpha *ABSUM1 + beta*C[(rowC+bdx)+ldc*colC];
	  }


	  if ((rowC < m) && (colC+bdy < n)) {
	    C[ rowC + ldc*(colC+bdy) ] = alpha *ABSUM2 + beta*C[rowC+ldc*(colC+bdy)];
	  }

	  if ((rowC+bdx < m) && (colC+bdy < n)) {
	    C[ (rowC+bdx) + ldc*(colC+bdy) ] = alpha *ABSUM3 + beta*C[(rowC+bdx)+ldc*(colC+bdy)];
	  }

	}
      }

      /* synchronize threads */
      __syncthreads();
      
      if ( tx == 0 && ty == 0 ) {
	s_count = atomicAdd(d_counter,1);
      } 
      __syncthreads();
      count = s_count;
      
    }

}



/* batched DGEMM kernel */
__global__ void dgemm_custom_simple_1block_batch_kernel_05( int transa, int transb,
		                                         int *mlist, int *nlist, int *klist,
                    			                 double alpha,
		                                         const double **Alist, int *ldalist,
                    			                 const double **Blist, int *ldblist,
		                                         double beta,
                    			                 double **Clist, int *ldclist,
							 int nbatch, 
							 int *d_counter)
{

    /* allocate shared memory */
  __shared__ double As[(blk_m+1)*blk_k];
  __shared__ double Bs[4][(blk_k+1)*blk_n];
  __shared__ int s_count;
  
  /* set thread index */
  int tx = threadIdx.x;
  int ty = threadIdx.y;

  int bdy = blockDim.y;
  int bdx = blockDim.x;

    /* local variables */ 
    int im, in, ik, i, id;
    int colA, colB, colC;
    int rowA, rowB, rowC;
    double ABSUM = 0.0;
 
    if ( tx == 0 && ty == 0 ) {
       s_count = atomicAdd(d_counter,1);
    } 

    __syncthreads();

    int count = s_count;

    if ( s_count >= nbatch ) return;

    /* get dimensions for current block (DGEMM) */
    int m = mlist[count];
    int n = nlist[count]; 
    int k = klist[count];
    int lda = ldalist[count];
    int ldb = ldblist[count];
    int ldc = ldclist[count];
    const double *A = Alist[count];
    const double *B = Blist[count];
    double *C = Clist[count];

    /*load L*/
    for ( ik=0; ik<k; ik+=bdy ) {
      for ( im=0; im<m; im+=bdx ) {
	if ((tx + im < m) && (ty + ik < k)) {
	  As[(tx+im)+(blk_m+1)*(ty+ik)] = __ldg(A+ (tx+im)+lda*(ty+ik) );
	}
      }
    }
    

    while ( s_count < nbatch ) {
    
      /* get dimensions for current block (DGEMM) */
      int m = mlist[count];
      int n = nlist[count]; 
      int k = klist[count];
      int lda = ldalist[count];
      int ldb = ldblist[count];
      int ldc = ldclist[count];
      const double *A = Alist[count];
      const double *B = Blist[count];
      double *C = Clist[count];

      __syncthreads();   
      
      /* early exit */
      if(!m || !n || !k) return;

      /* load B */
      for ( ik=0; ik<n; ik+=bdy ) {
	for ( im=0; im<k; im+=bdx ) {
	  if ((tx + im < k) && (ty + ik < n)) {
	    Bs[0][(tx+im)+(blk_k+1)*(ty+ik)] = __ldg(B+(tx+im)+ldb*(ty+ik) );
	  }
	}
      }

      /* synchronize threads */
      __syncthreads();
    
      for ( in=0; in<n; in+=bdy ) {
	for ( im=0; im<k; im+=bdx ) {

	  rowC = tx+im;
	  colC = ty+in;

	  if ((rowC < m) && (colC < n)) {
	    
	    ABSUM = 0.0;
	
	    /* compute A*B */
	    for(i = 0; i < k; i++){
	      ABSUM += As[rowC+(blk_m+1)*i]*Bs[0][i+(blk_k+1)*colC];
	    }
	
	    C[ rowC + ldc*colC ] = alpha *ABSUM + beta*C[rowC+ldc*colC];
	  }
	}
      }

      /* synchronize threads */
      __syncthreads();
      
      if ( tx == 0 && ty == 0 ) {
	s_count = atomicAdd(d_counter,1);
      } 
      __syncthreads();
      count = s_count;
      
    }

}



/* batched DGEMM kernel */
__global__ void dgemm_custom_simple_1block_batch_kernel_03( int transa, int transb,
		                                         int *mlist, int *nlist, int *klist,
                    			                 double alpha,
		                                         const double **Alist, int *ldalist,
                    			                 const double **Blist, int *ldblist,
		                                         double beta,
                    			                 double **Clist, int *ldclist,
							 int nbatch, 
							 int *d_counter)
{

    /* allocate shared memory */
  __shared__ double As[blk_m*(blk_k+1)];
  __shared__ double B0s[blk_k*(blk_n+1)];
  __shared__ double B1s[blk_k*(blk_n+1)];
  __shared__ double B2s[blk_k*(blk_n+1)];
  __shared__ double B3s[blk_k*(blk_n+1)];
    __shared__ int s_count0;
    __shared__ int s_count1;
    __shared__ int s_count2;
    __shared__ int s_count3;

    double *Bp;

    /* set thread index */
    const  int tx = threadIdx.x;
    const  int ty = threadIdx.y;

    /* local variables */ 
    int im, in, ik, i, id;
    int colA, colB, colC;
    int rowA, rowB, rowC;
    double ABSUM00 = 0.0;
    double ABSUM01 = 0.0;
    double ABSUM10 = 0.0;
    double ABSUM11 = 0.0;
    int rowC1, colC1;
 
    if ( tx ==0 && ty == 0 ) printf ("entering \n" );

    if ( tx == 0 && ty == 0 ) {
       s_count0 = atomicAdd(d_counter,1);
    } 

    if ( tx == 8 && ty == 0 ) {
       s_count1 = atomicAdd(d_counter,1);
    } 

    if ( tx == 0 && ty == 8 ) {
       s_count2 = atomicAdd(d_counter,1);
    } 

    if ( tx == 8 && ty == 8 ) {
       s_count3 = atomicAdd(d_counter,1);
    } 

    __syncthreads();

    int count;
    if ( tx/8 == 0 && ty/8 == 0 )  {
      count = s_count0;
      Bp = B0s;
    }
    if ( tx/8 == 1 && ty/8 == 0 )  {
      count = s_count1;
      Bp = B1s;
    }
    if ( tx/8 == 0 && ty/8 == 1 )  {
      count = s_count2;
      Bp = B2s;
    }
    if ( tx/8 == 1 && ty/8 == 1 )  {
      count = s_count3;
      Bp = B3s;
    }

    if ( tx ==0 && ty == 0 ) printf ("first count = %d \n", count );

    if ( count >= nbatch ) return;

    /* get dimensions for current block (DGEMM) */
    int m = mlist[count];
    int n = nlist[count]; 
    int k = klist[count];
    int lda = ldalist[count];
    int ldb = ldblist[count];
    int ldc = ldclist[count];
    const double *A = Alist[count];
    const double *B = Blist[count];
    double *C = Clist[count];
    double A0, A1, B0, B1;

    /*load L*/
    for ( ik=0; ik<k; ik+=blockDim.y ) {
      for ( im=0; im<m; im+=blockDim.x ) {
	if ((tx + im < m) && (ty + ik < k)) {
	  As[(tx+im)+(ty+ik)*lda] = __ldg(A+ (tx+im) + lda*(ty+ik) );
	}
      }
    }

    if ( tx ==0 && ty == 0 ) printf ("l loaded \n" );

    while ( count < nbatch ) {
    
      /* get dimensions for current block (DGEMM) */
      int m = mlist[count];
      int n = nlist[count]; 
      int k = klist[count];
      int lda = ldalist[count];
      int ldb = ldblist[count];
      int ldc = ldclist[count];
      const double *A = Alist[count];
      const double *B = Blist[count];
      double *C = Clist[count];

      __syncthreads();   
      
      /* early exit */
      if(!m || !n || !k) return;

      /* load B - just using 8x8 threads*/
      
      for ( ik=0; ik<n; ik+=blockDim.y/2 ) {
	for ( im=0; im<k; im+=blockDim.x/2 ) {
	  if ((tx%8 + im < k) && (ty%8 + ik < n)) {
	    Bp[tx%8+im + (ty%8+ik)*ldb] = __ldg(B+(tx%8+im) + ldb*(ty%8+ik) );
	  }
	}
      }
      

      /* synchronize threads */
      __syncthreads();

      for ( in=0; in<n; in+=blockDim.y ) {
	for ( im=0; im<k; im+=blockDim.x ) {
	  rowC = tx%8+im;
	  colC = ty%8+in;
	  rowC1 = rowC + blockDim.x/2;
	  colC1 = colC + blockDim.y/2;

	  ABSUM00 = 0.0;
	  ABSUM01 = 0.0;
	  ABSUM10 = 0.0;
	  ABSUM11 = 0.0;
	  
	  for ( i=0; i<k; i++ ) {
	    
	    // read from shared //
	    if ( rowC < m ) A0 = As[rowC+lda*i];
	    if ( rowC1 < m ) A1 = As[rowC1+lda*i];
	    if ( colC < n ) B0 = Bp[i+ldb*colC];
	    if ( colC1 < n ) B1 = Bp[i+ldb*colC1];
	    
	    if ((rowC < m) && (colC < n)) {
	      ABSUM00 += A0*B0;
	    }
	    if ((rowC1 < m) && (colC < n)) {
	      ABSUM10 += A1*B0;
	    }
	    if ((rowC < m) && (colC1 < n)) {
	      ABSUM01 += A0*B1;
	    }
	    if ((rowC1 < m) && (colC1 < n)) {
	      ABSUM11 += A1*B1;
	    }
	    
	  }

	  if ((rowC < m) && (colC < n)) {
	    C[rowC + ldc*colC] = alpha * ABSUM00 + beta * C[rowC+ldc*colC];
	  }
	  
	  if ((rowC1 < m) && (colC < n)) {
	    C[rowC1 + ldc*colC] = alpha * ABSUM10 + beta * C[rowC1+ldc*colC];
	  }
	  
	  if ((rowC < m) && (colC1 < n)) {
	    C[rowC + ldc*colC1] = alpha * ABSUM01 + beta * C[rowC+ldc*colC1];
	  }
	  
	  if ((rowC1 < m) && (colC1 < n)) {
	    C[rowC1 + ldc*colC1] = alpha * ABSUM11 + beta * C[rowC1+ldc*colC1];
	  }	  

	}
      }

      /* synchronize threads */
      __syncthreads();

      /*
      if ( tx == 0 && ty == 0 ) {
	s_count = atomicAdd(d_counter,1);
      } 
      __syncthreads();
      count = s_count;
      */
      if ( tx == 0 && ty == 0 ) {
	s_count0 = atomicAdd(d_counter,1);
	printf ("scount0 = %d, nbatch = %d \n", s_count0, nbatch);
      } 

      if ( tx == 8 && ty == 0 ) {
	s_count1 = atomicAdd(d_counter,1);
      } 

      if ( tx == 0 && ty == 8 ) {
	s_count2 = atomicAdd(d_counter,1);
      } 
      
      if ( tx == 8 && ty == 8 ) {
	s_count3 = atomicAdd(d_counter,1);
      } 
      
      __syncthreads();

      int count;
      if ( tx/8 == 0 && ty/8 == 0 )  {
	count = s_count0;
	Bp = B0s;
      }
      if ( tx/8 == 1 && ty/8 == 0 )  {
	count = s_count1;
	Bp = B1s;
      }
      if ( tx/8 == 0 && ty/8 == 1 )  {
	count = s_count2;
	Bp = B2s;
      }
      if ( tx/8 == 1 && ty/8 == 1 )  {
	count = s_count3;
	Bp = B3s;
      }
      
    }

    return;

}


/* batched DGEMM kernel */
__global__ void dgemm_custom_simple_1block_batch_kernel_02_working_200GF( int transa, int transb,
		                                         int *mlist, int *nlist, int *klist,
                    			                 double alpha,
		                                         const double **Alist, int *ldalist,
                    			                 const double **Blist, int *ldblist,
		                                         double beta,
                    			                 double **Clist, int *ldclist,
							 int nbatch, 
							 int *d_counter)
{

    /* allocate shared memory */
    __shared__ double As[blk_m][blk_k+1];
    __shared__ double Bs[blk_k][blk_n+1];
    __shared__ int s_count;

    /* set thread index */
    const  int tx = threadIdx.x;
    const  int ty = threadIdx.y;

    /* local variables */ 
    int im, in, ik, i, id;
    int colA, colB, colC;
    int rowA, rowB, rowC;
    double ABSUM = 0.0;
 
    if ( tx == 0 && ty == 0 ) {
       s_count = atomicAdd(d_counter,1);
    } 

    __syncthreads();

    int count = s_count;

    if ( s_count >= nbatch ) return;

    /* get dimensions for current block (DGEMM) */
    int m = mlist[count];
    int n = nlist[count]; 
    int k = klist[count];
    int lda = ldalist[count];
    int ldb = ldblist[count];
    int ldc = ldclist[count];
    const double *A = Alist[count];
    const double *B = Blist[count];
    double *C = Clist[count];

    /*load L*/
    for ( ik=0; ik<k; ik+=blockDim.y ) {
      for ( im=0; im<m; im+=blockDim.x ) {
	if ((tx + im*blockDim.x < m) && (ty + ik*blockDim.y < k)) {
	  As[tx+im*blockDim.x][ty+ik*blockDim.y] = __ldg(A+ (tx+im*blockDim.x) + lda*(ty+ik*blockDim.y) );
	}
      }
    }
    

    while ( s_count < nbatch ) {
    
      /* get dimensions for current block (DGEMM) */
      int m = mlist[count];
      int n = nlist[count]; 
      int k = klist[count];
      int lda = ldalist[count];
      int ldb = ldblist[count];
      int ldc = ldclist[count];
      const double *A = Alist[count];
      const double *B = Blist[count];
      double *C = Clist[count];

      __syncthreads();   
      
      /* early exit */
      if(!m || !n || !k) return;

      /* load B */
      for ( ik=0; ik<n; ik+=blockDim.y ) {
	for ( im=0; im<k; im+=blockDim.x ) {
	  if ((tx + im*blockDim.x < k) && (ty + ik*blockDim.y < n)) {
	    Bs[tx+im*blockDim.x][ty+ik*blockDim.y] = __ldg(B+(tx+im*blockDim.x) + ldb*(ty+ik*blockDim.y) );
	  }
	}
      }

      /* synchronize threads */
      __syncthreads();
    
      for ( in=0; in<n; in+=blockDim.y ) {
	for ( im=0; im<k; im+=blockDim.x ) {
	  rowC = tx+im*blockDim.x;
	  colC = ty+in*blockDim.y;
	  if ((rowC < m) && (colC < n)) {
	    
	    ABSUM = 0.0;
	
	    /* compute A*B */
	    for(i = 0; i < k; i++){
	      ABSUM += As[rowC][i]*Bs[i][colC];
	    }
	
	    C[ rowC + ldc*colC ] = alpha *ABSUM + beta*C[rowC+ldc*colC];
	  }
	}
      }

      /* synchronize threads */
      __syncthreads();
      
      if ( tx == 0 && ty == 0 ) {
	s_count = atomicAdd(d_counter,1);
      } 
      __syncthreads();
      count = s_count;
      
    }

}


/* batched DGEMM kernel */
__global__ void dgemm_custom_simple_1block_batch_kernel_01( int transa, int transb,
		                                         int *mlist, int *nlist, int *klist,
                    			                 double alpha,
		                                         const double **Alist, int *ldalist,
                    			                 const double **Blist, int *ldblist,
		                                         double beta,
                    			                 double **Clist, int *ldclist,
							 int nbatch, 
							 int *d_counter)
{

    /* allocate shared memory */
    __shared__ double As[blk_m][blk_k+1];
    __shared__ double Bs[blk_k][blk_n+1];
    __shared__ int s_count;

    /* set thread index */
    const  int tx = threadIdx.x;
    const  int ty = threadIdx.y;

    /* local variables */ 
    int im, in, ik, i, id;
    int colA, colB;
    int rowA, rowB;
    double ABSUM = 0.0;
 
    if ( tx == 0 && ty == 0 ) {
       s_count = atomicAdd(d_counter,1);
    } 

    __syncthreads();

    int count = s_count;

    if ( s_count >= nbatch ) return;

    /* get dimensions for current block (DGEMM) */
    int m = mlist[count];
    int n = nlist[count]; 
    int k = klist[count];
    int lda = ldalist[count];
    int ldb = ldblist[count];
    int ldc = ldclist[count];
    const double *A = Alist[count];
    const double *B = Blist[count];
    double *C = Clist[count];


    /*load L*/

    rowA = tx;
    colA = ty;

    if((rowA < m) && (colA < k)){     
      As[tx][ty] = __ldg(A+(rowA + lda*colA));                 
    }
    

    while ( s_count < nbatch ) {
    
      /* get dimensions for current block (DGEMM) */
      int m = mlist[count];
      int n = nlist[count]; 
      int k = klist[count];
      int lda = ldalist[count];
      int ldb = ldblist[count];
      int ldc = ldclist[count];
      const double *A = Alist[count];
      const double *B = Blist[count];
      double *C = Clist[count];

      __syncthreads();   
      
      /* early exit */
      if(!m || !n || !k) return;

      /* load B */
      if ((tx < k) && (ty < n)) {
	Bs[tx][ty] = __ldg(B+(tx + ldb*ty));          
      }

      /* synchronize threads */
      __syncthreads();
    
      if ( tx < m && ty < n ) {

	ABSUM = 0.0;
	
	/* compute A*B */
	for(i = 0; i < k; i++){
	  ABSUM += As[tx][i]*Bs[i][ty];
	}
	
	C[ tx + ldc*ty ] = alpha *ABSUM + beta*C[tx+ldc*ty];
      }

      /* synchronize threads */
      __syncthreads();
      
      if ( tx == 0 && ty == 0 ) {
	s_count = atomicAdd(d_counter,1);
      } 
      __syncthreads();
      count = s_count;
      
    }

}


