/******************************************************************************
 *
 * Routines for doing timing.
 *
 *****************************************************************************/

#include "transport_headers.h"

/*-------------------------------------------------------
 * Timing macros
 *-------------------------------------------------------*/

#define StartTiming() \
  TimingWallCount -= MPI_Wtime();       \
  TimingCPUCount -= 0.0
/*time_getCPUSeconds()*/

#define StopTiming() \
  TimingWallCount += MPI_Wtime();       \
  TimingCPUCount += 0.0
/*time_getCPUSeconds()*/

#define global_timing_ref(field) global_timing->field

TimingType *global_timing = NULL;

/*--------------------------------------------------------------------------
 * InitializeTiming
 *--------------------------------------------------------------------------*/

int
InitializeTiming( const char *name )
{
  int timing_index;

  double  *old_wall_time;
  double  *old_cpu_time;
  double  *old_flops;
  char   **old_name;
  int     *old_state;
  int     *old_num_regs;

  int new_name;
  int i;

  /*-------------------------------------------------------
   * Allocate global TimingType structure if needed
   *-------------------------------------------------------*/

  if(global_timing == NULL){
    CNEW(global_timing, 1, TimingType *);
  }

  /*-------------------------------------------------------
   * Check to see if name has already been registered
   *-------------------------------------------------------*/

  new_name = 1;
  for(i = 0; i < (global_timing_ref(size)); i++){
    if(TimingNumRegs(i) > 0){
      if(strcmp(name, TimingName(i)) == 0){
        new_name = 0;
        timing_index = i;
        TimingNumRegs(timing_index)++;
        break;
      }
    }
  }

  if(new_name){
    for(i = 0; i < global_timing_ref(size); i++){
      if(TimingNumRegs(i) == 0){
        break;
      }
    }
    timing_index = i;
  }

  /*-------------------------------------------------------
   * Register the new timing name
   *-------------------------------------------------------*/

  if(new_name){
    if(timing_index == (global_timing_ref(size))){
      old_wall_time = (global_timing_ref(wall_time));
      old_cpu_time  = (global_timing_ref(cpu_time));
      old_flops     = (global_timing_ref(flops));
      old_name      = (global_timing_ref(name));
      old_state     = (global_timing_ref(state));
      old_num_regs  = (global_timing_ref(num_regs));

      CNEW(global_timing_ref(wall_time), timing_index+1, double *);
      CNEW(global_timing_ref(cpu_time), timing_index+1, double *);
      CNEW(global_timing_ref(flops), timing_index+1, double *);
      CNEW(global_timing_ref(name), timing_index+1, char **);
      CNEW(global_timing_ref(state), timing_index+1, int *);
      CNEW(global_timing_ref(num_regs), timing_index+1, int *);
      (global_timing_ref(size))++;

      for(i = 0; i < timing_index; i++){
        TimingWallTime(i) = old_wall_time[i];
        TimingCPUTime(i)  = old_cpu_time[i];
        TimingFLOPS(i)    = old_flops[i];
        TimingName(i)     = old_name[i];
        TimingState(i)    = old_state[i];
        TimingNumRegs(i)  = old_num_regs[i];
      }

      FREE(old_wall_time);
      FREE(old_cpu_time);
      FREE(old_flops);
      FREE(old_name);
      FREE(old_state);
      FREE(old_num_regs);
    }

    CNEW(TimingName(timing_index), 80, char *);
    strncpy(TimingName(timing_index), name, 79);
    TimingState(timing_index)   = 0;
    TimingNumRegs(timing_index) = 1;
    (global_timing_ref(num_names))++;
  }

  return( timing_index);
}

/*--------------------------------------------------------------------------
 * FinalizeTiming
 *--------------------------------------------------------------------------*/

int
FinalizeTiming( int timing_index )
{
  int ierr = 0;

  if(global_timing == NULL){
    return( ierr);
  }

  if(timing_index < (global_timing_ref(size))){
    if(TimingNumRegs(timing_index) > 0){
      TimingNumRegs(timing_index)--;
    }

    if(TimingNumRegs(timing_index) == 0){
      FREE(TimingName(timing_index));
      (global_timing_ref(num_names))--;
    }
  }

  if((global_timing->num_names) == 0){
    FREE(global_timing_ref(wall_time));
    FREE(global_timing_ref(cpu_time));
    FREE(global_timing_ref(flops));
    FREE(global_timing_ref(name));
    FREE(global_timing_ref(state));
    FREE(global_timing_ref(num_regs));

    FREE(global_timing);
    global_timing = NULL;
  }

  return( ierr);
}

/*--------------------------------------------------------------------------
 * IncFLOPCount
 *--------------------------------------------------------------------------*/

int
IncFLOPCount( int inc )
{
  int ierr = 0;

  if(global_timing == NULL){
    return( ierr);
  }

  TimingFLOPCount += (double) (inc);

  return( ierr);
}

/*--------------------------------------------------------------------------
 * BeginTiming
 *--------------------------------------------------------------------------*/

int
BeginTiming( int timing_index )
{
  int ierr = 0;

  if(global_timing == NULL){
    return( ierr);
  }

  if(TimingState(timing_index) == 0){
    StopTiming();
    TimingWallTime(timing_index) -= TimingWallCount;
    TimingCPUTime(timing_index)  -= TimingCPUCount;
    TimingFLOPS(timing_index)    -= TimingFLOPCount;

    StartTiming();
  }
  TimingState(timing_index)++;

  return( ierr);
}

/*--------------------------------------------------------------------------
 * EndTiming
 *--------------------------------------------------------------------------*/

int
EndTiming( int timing_index )
{
  int ierr = 0;

  if(global_timing == NULL){
    return( ierr);
  }

  TimingState(timing_index)--;
  if(TimingState(timing_index) == 0){
    StopTiming();
    TimingWallTime(timing_index) += TimingWallCount;
    TimingCPUTime(timing_index)  += TimingCPUCount;
    TimingFLOPS(timing_index)    += TimingFLOPCount;
    StartTiming();
  }

  return( ierr);
}

/*--------------------------------------------------------------------------
 * ClearTiming
 *--------------------------------------------------------------------------*/

int
ClearTiming( )
{
  int ierr = 0;
  int i;

  if(global_timing == NULL){
    return( ierr);
  }

  for(i = 0; i < (global_timing_ref(size)); i++){
    TimingWallTime(i) = 0.0;
    TimingCPUTime(i)  = 0.0;
    TimingFLOPS(i)    = 0.0;
  }

  return( ierr);
}

/*--------------------------------------------------------------------------
 * PrintTiming
 *--------------------------------------------------------------------------*/

int
PrintTiming( const char     *heading,
             MPI_Comm comm  )
{
  int ierr = 0;

  double local_wall_time;
  double local_cpu_time;
  double wall_time;
  double cpu_time;
  double wall_mflops;
  double cpu_mflops;

  int i;
  int myrank;

  if(global_timing == NULL){
    return( ierr);
  }

  MPI_Comm_rank(comm, &myrank );

  /* print heading */
  if(myrank == 0){
    printf("\n");
    printf("============================================================\n");
    printf("%s:\n", heading);
    printf("============================================================\n");
  }

  for(i = 0; i < (global_timing->size); i++){
    if(TimingNumRegs(i) > 0){
      local_wall_time = TimingWallTime(i);
      local_cpu_time  = TimingCPUTime(i);
      MPI_Allreduce(&local_wall_time, &wall_time, 1,
                    MPI_DOUBLE, MPI_MAX, comm);
      MPI_Allreduce(&local_cpu_time, &cpu_time, 1,
                    MPI_DOUBLE, MPI_MAX, comm);

      if(myrank == 0){
        printf("%s:\n", TimingName(i));

        /* print wall clock info */
        printf("  wall clock time = %f seconds\n", wall_time);
        if(wall_time){
          wall_mflops = TimingFLOPS(i) / wall_time / 1.0E6;
        }
        else {
          wall_mflops = 0.0;
        }
        /*printf("  wall MFLOPS     = %f\n", wall_mflops);*/

        /* print CPU clock info */
        /*
            printf("  cpu clock time  = %f seconds\n", cpu_time);
            if(cpu_time){
              cpu_mflops = TimingFLOPS(i) / cpu_time / 1.0E6;
            }
            else {
              cpu_mflops = 0.0;
          }*/
        /*printf("  cpu MFLOPS      = %f\n\n", cpu_mflops);*/
      }
    }
  }

  return( ierr);
}
