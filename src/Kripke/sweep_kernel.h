/******************************************************************************
 *
 * Header for sweep_kernel
 *
 *****************************************************************************/

#ifndef _SWEEP_KERNEL_HEADER
#define _SWEEP_KERNEL_HEADER

#include "transport_headers.h"

/*--------------------------------------------------------------------------
 * Structures
 *--------------------------------------------------------------------------*/

typedef struct {
  MPI_Comm comm;
  Sweep_Solver_Data    *sweep_solver_data;
} Boltzmann_Solver;

#endif
