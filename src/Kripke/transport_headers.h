/*--------------------------------------------------------------------------
 * Header files for 1-D transport code.
 *--------------------------------------------------------------------------*/

#ifndef transport_header
#define transport_header

/*--------------------------------------------------------------------------
 * Include standard C libraries
 *--------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

/*--------------------------------------------------------------------------
 * Include data structures for sweep kernel code
 *--------------------------------------------------------------------------*/

#include "bc_data.h"
#include "data_vector.h"
#include "grid.h"
#include "input_variables.h"
#include "sigma_tot.h"
#include "sweep_solver.h"
#include "sweep_kernel.h"
#include "user_data.h"
#include "timing.h"

/*--------------------------------------------------------------------------
 * Include all function prototypes.
 *--------------------------------------------------------------------------*/

#include "transport_protos.h"
#include "comm.h"
#include<stdio.h>
#include<stdlib.h>

#define NEW(ptr, nitems, \
            cast) if(((ptr) = \
                        (cast) malloc((unsigned) ((nitems) * \
                                                  sizeof(*(ptr))))) == \
                     NULL){fprintf( \
                             stderr, \
                             "ERROR: malloc failed, file %s, line %d\n", \
                             __FILE__, \
                             __LINE__); error_exit(1); }

#define CNEW(ptr, nitems, \
             cast) if(((ptr) = \
                         (cast) calloc((nitems), \
                                       (unsigned) (sizeof(*(ptr))))) == \
                      NULL){fprintf( \
                              stderr, \
                              "ERROR: malloc failed, file %s, line %d\n", \
                              __FILE__, \
                              __LINE__); \
                            error_exit( \
                              1); }

#define FREE(ptr)  {void *ptr_temp; ptr_temp = ptr; free( (ptr_temp) ); }

#ifndef MIN
#define MIN(A, B) ((A) < (B) ? (A) : (B))
#endif

#ifndef MAX
#define MAX(A, B) ((A) > (B) ? (A) : (B))
#endif

#ifndef ABS
#define ABS(A)    ((A < 0) ? -(A) : (A))
#endif

#endif
