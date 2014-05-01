/******************************************************************************
 *
 * Header file for doing timing
 *
 *****************************************************************************/

#ifndef TIMING_HEADER
#define TIMING_HEADER

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#ifdef __cplusplus
extern "C" {
#endif

/*--------------------------------------------------------------------------
 * Prototypes for low-level timing routines
 *--------------------------------------------------------------------------*/

/* timer.c */
/*double time_getWallclockSeconds( void );
double time_getCPUSeconds( void );
double time_get_wallclock_seconds_( void );
double time_get_cpu_seconds_( void );*/

/*-------------------------------------------------------
 * Global timing structure
 *-------------------------------------------------------*/

typedef struct {
  double  *wall_time;
  double  *cpu_time;
  double  *flops;
  char   **name;
  int     *state;      /* boolean flag to allow for recursive timing */
  int     *num_regs;   /* count of how many times a name is registered */

  int num_names;
  int size;

  double wall_count;
  double CPU_count;
  double FLOP_count;

} TimingType;

/*-------------------------------------------------------
 * Accessor functions
 *-------------------------------------------------------*/

#define TimingWallTime(i) (global_timing->wall_time[(i)])
#define TimingCPUTime(i)  (global_timing->cpu_time[(i)])
#define TimingFLOPS(i)    (global_timing->flops[(i)])
#define TimingName(i)     (global_timing->name[(i)])
#define TimingState(i)    (global_timing->state[(i)])
#define TimingNumRegs(i)  (global_timing->num_regs[(i)])
#define TimingWallCount   (global_timing->wall_count)
#define TimingCPUCount    (global_timing->CPU_count)
#define TimingFLOPCount   (global_timing->FLOP_count)

/*-------------------------------------------------------
 * Prototypes
 *-------------------------------------------------------*/

/* timing.c */
int InitializeTiming( const char *name );
int FinalizeTiming( int timing_index );
int IncFLOPCount( int inc );
int BeginTiming( int timing_index );
int EndTiming( int timing_index );
int ClearTiming( void );
int PrintTiming( const char *heading, MPI_Comm comm );

enum TimingTypes {SWEEP=0, SCATTERING, FISSION, LTIMES, LPLUSTIMES,
                  CREATEPQRGRID};

#ifdef __cplusplus
}
#endif

#endif
