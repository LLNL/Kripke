/*

   This is the header file containing the function prototypes
   for the communication layer implemented atop MPI.  This header
   needs to be included in any file containing functions that call
   functions in the communication layer.  In effect, the function
   prototypes in {\bf comm.h} define the communication layer.  That
   is, if the layer needs to be ported to a new machine or
   communication system, all of the functions declared in this
   prototype must be implemented.

*/

#include <Kripke/Comm.h>
#include <fcntl.h>
#include <unistd.h>
#include <mpi.h>
#include <vector>
#include <stdio.h>


/**
   Initializes storage for messages (and associated status arrays)
   posted in one of the 6 signed coordinate directions.  The signed coordinate
   directions are enumerated by a key as follows:

            key      direction
             0       positive x
             1       negative x
             2       positive y
             3       negative y
             4       positive z
             5       negative z

   The parameters ``len'' and ``nm'' are each pointers to integer arrays of
   length 6.  For key = 0, ..., 5, len[key] is the length (number of doubles)
   of each message to be posted in the direction associated with key, and
   nm[key] is the number of such messages being posted in the associated
   direction.
*/
Comm::Comm(int * len, int *  nm){
  int i, j, k;
  int size;

  size = 0;
  buf_total = 0;
  for(i=0; i<6; i++){
    which_len[i] = len[i];
    which_num[i] = nm[i];
    size += len[i] * nm[i];   /* space for which direction bufs*/
    buf_total += nm[i];
  }

  int mpi_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  if(mpi_rank == 0){
    printf("Comm: allocating %d buffers (%lf megabytes)\n", buf_total, (double)(size*sizeof(double))/(1024.0*1024.0));
  }

  buf_which[0] = 0;
  for(i=1; i<6; i++){
    buf_which[i] = buf_which[i-1] + which_num[i-1];
  }

  bptr_sv.resize(size);
  double *bptr = &bptr_sv[0];

  buf_pool.resize(buf_total);

  for(i = 0; i < 6; i++){    /* for each "which" set */
    k = which_len[i];       /* msg size in this set */
    for(j = 0; j < which_num[i]; j++){  /* for each msg in the set */
      buf_pool[buf_which[i] + j] = bptr;     /* ptr to next message */
      bptr = bptr + k;    /* bump over this message */
    }
  }

  buf_s_req.resize(buf_total);
  buf_r_req.resize(buf_total);
  buf_status.resize(buf_total);

  buf_reset();
}

/**
   R_recv_dir() posts receives in the signed direction ``which''.  The number
   of receives posted is nm[which], where nm[] is the array passed to
   the constructor of Comm.  The argument ``r_member''
   is the r coordinate of the node from which the nm[which] messages will
   be sent.  R_recv_dir is non-blocking, i.e., it returns after posting
   the receives without waiting for any messages to arrive.
*/
void Comm::R_recv_dir(int which, int r_member){
  int i, j;

  i = buf_which[which];     /* get to right buffer pool */
  for(j = 0; j < which_num[which]; j++){
    MPI_Irecv(buf_pool[i+j], which_len[which], MPI_DOUBLE, r_member,
              0, MPI_COMM_WORLD, &buf_r_req[i+j]);
  }

  buf_rec[which]+= 1;  /*record recv in process: += is debug error guard */
}

/**
   R_recv_test() checks to see if any of the posted receives in the
   signed direction ``which'' has been satisfied (i.e, a message has
   arrived).  If no received message is pending, the function returns 0.
   If a message has arrived, its address is assigned to *msg and the
   function returns 1.
*/
int Comm::R_recv_test(int which,  double **msg){


  int i, j;
  int done;      /* index of which finished, if any */
  int flag;      /* 0 if none - 1 if someone finished */
  MPI_Status status;    /* status of one that finished, if any */

  i = buf_which[which];
  j = buf_rec[which];

  if(j > 0){     /* non-blocking rcvS posted in this direction - any done?*/
    int k, reqs_remain=0;

    /* the following loop tests to see if there are any pending requests
       before calling MPI_Testany.  The 1.2a version of the Edinburgh
       T3D MPI port of MPI_Testany returned a true flag with a -1
       index value if a list of null requests was passed to it
       (ie. a bug)
       */
    for(k=0; k<which_num[which]; k++){
      if(buf_r_req[i+k] != MPI_REQUEST_NULL){
        reqs_remain++;
      }
    }

    if(reqs_remain > 0){
      MPI_Testany(which_num[which], &buf_r_req[i], &done, &flag,
                  &status);
    }
    else {
      flag = 0;
    }

    if(flag > 0){
      *msg = buf_pool[i + done];
    }
    return( flag);      /*  note - msg is returned only with flag = 1 */
  }
  /* this direction is boundary - get  next empty buffer  == new*/
  *msg = buf_pool[i-j];        /* remember j is <=0 , so this is an add */
  buf_rec[which] -= 1;         /* bump to bigger negative for next call */
  return( 1);
}


/**
  post non-blocking send & remember handle to check completion later
  send-handles are all mixed together in one array
  Called from a node (p,q,r), R_send() sends the message of ``length''
  doubles pointed at by ``msg'' to node (p,q,r_member).  The send is
  non-blocking, i.e., it returns immediately after submitting the message
  to the underlying communication system.
*/
void Comm::R_send(double * msg, int r_member, int length){
  if(r_member >= 0){
    MPI_Isend(msg, length, MPI_DOUBLE, r_member, 0, MPI_COMM_WORLD,
              &buf_s_req[send_cnt]);
    send_cnt++; /* bump for next send */
  }
}


/**
   R_wait_send() waits for all posted sends to complete.
*/
void Comm::R_wait_send (){
  MPI_Waitall(send_cnt, &buf_s_req[0], &buf_status[0]);

  buf_reset();    /* all done - reset all the buffer info */
}



/**
 * initialize all the things that change each iteration
 */
void Comm::buf_reset(void){
  send_cnt = 0;

  int i;
  for(i = 0; i < 6; i++){
    buf_rec[i] = 0;
  }
  for(i = 0; i < buf_total; i++){
    buf_r_req[i] = MPI_REQUEST_NULL;
    buf_s_req[i] = MPI_REQUEST_NULL;
  }
}




/*========================= MISCELLANEOUS ======================*/

/**
   error_exit() aborts the concurrent execution, passing the flag value
   to lower-level (system-dependent) abort routines, which may ignore it.
*/
void error_exit(int flag)
{
  fflush(NULL);
  MPI_Abort(MPI_COMM_WORLD, flag);
}
