/******************************************************************************
 *
 * Header for sweep_solver
 *
 *****************************************************************************/

#ifndef _SWEEP_SOLVER_HEADER
#define _SWEEP_SOLVER_HEADER

#include "transport_headers.h"

/*--------------------------------------------------------------------------
 * Structures
 *--------------------------------------------------------------------------*/

typedef struct {
  int                   *swept;
  int                   *swept_min;
  int                   *swept_max;

  MPI_Comm comm;

} Sweep_Solver_Data;

#endif
