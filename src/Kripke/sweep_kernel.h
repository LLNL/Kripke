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

struct Boltzmann_Solver {
  MPI_Comm comm;
  Sweep_Solver_Data    *sweep_solver_data;
};

#endif
