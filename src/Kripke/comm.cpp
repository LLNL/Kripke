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

#include <fcntl.h>
#include <unistd.h>
#include "transport_headers.h"
#include "comm.h"


/* declarations for MPI version of the code */
MPI_Comm R_group;

/* Saved values of R and r */
static int R0;
static int r0;

/*  -------------------------- message buffer declarations -------------- */
/* message buffer information, used in send-receive calls */
/* the length of several arrays is not known until runtime, */
/* so will be allocated via malloc in the buf_init routine   */

/* these items are initialized in buf_init and don't change */
int buf_total;     /* count of total number of messages */
int which_num[6];  /* max number of message in each direction */
int which_len[6];  /* size of messages in each direction */
int buf_which[6];  /* array to say where in pool each direction starts */
double ** buf_pool; /* array of  pointers to each buffer in pool */
double * bptr_sv;  /* pointer to memory for all buffers */

/* these items are reinitialized every iteration through the main code */
int send_cnt;     /* index into buf_s_req for next send-handle */
int buf_rec[6];   /* array to track whether which-direction is boundary */
                  /* if value is 1 - receives have been posted in direction */
                  /* if values is > 1 error: rcv's posted more than once*/
                  /* if value is <= 0 direction may be used for boundry bufs*/
                  /* amount less than  0 tracks which to allocate next  */
MPI_Request *buf_r_req; /* array to hold async recv message handles */
                        /* entries 1-1 with buf_pool, initialed to null-handle
                          */
MPI_Request *buf_s_req; /* array to hold async send message handles */
                        /* entries made in order send's posted (see send_cnt)
                          */
MPI_Status *buf_status; /* not used except for debug and MPI returns it */

void buf_reset();  /* grump - decl just to avoid an error message */



/*================= POINT TO POINT MESSAGE PASSING ==================*/

/* note previous versions had a new and free call for use with boundaries
     */
/* in this mpi version --- the "new" function is handled via R_recv_test */
/* and the "free" function is handled via reinitialization after R_send_wait
  */

void
R_recv_dir(int which,
           int r_member)
/* post message receives - all receives in specific direction posted at once
  */

{
  int i, j;

  i = buf_which[which];     /* get to right buffer pool */
  for(j = 0; j < which_num[which]; j++){
    MPI_Irecv(buf_pool[i+j], which_len[which], MPI_DOUBLE, r_member,
              0, R_group, &buf_r_req[i+j]);
  }

  buf_rec[which]+= 1;  /*record recv in process: += is debug error guard */
}

int
R_recv_test(int which,
            double **msg)
{  /* test if any message buffer is complete and ready to use */
   /* return 0 if not; return 1 and *msg if yes */

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

void
R_send(double * msg,
       int r_member,
       int length)
{
  /* post non-blocking send & remember handle to check completion later */
  /* send-handles are all mixed together in one array */

  if(r_member >= 0){
    MPI_Isend(msg, length, MPI_DOUBLE, r_member, 0, R_group,
              &buf_s_req[send_cnt]);
    send_cnt++; /* bump for next send */
  }
}

void
R_wait_send ()
{
  /* come here to check that all sends have completed so that the */
  /* buffers are free for reuse on next iteration */

  MPI_Waitall(send_cnt, buf_s_req, buf_status);

  buf_reset();    /* all done - reset all the buffer info */
}


/*========================= MISCELLANEOUS ======================*/

void
error_exit(int flag)
{
  fflush(NULL);
  MPI_Abort(MPI_COMM_WORLD, flag);
}

MPI_Comm GetRGroup() {
  return(R_group);
}

int GetRrank() {
  return(r0);
}

/*===================INITIALIZATION AND TERMINATION =============*/

void
create_R_grid(int R)
{
  int myrank, p, q, r, size;

  /* Load R into saved version. */
  R0 = R;

  /* Gte rank */
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  r = myrank;

  /* Load r into saved version. */
  r0 = r;

  /*------------------------R Group---------------------------*/
  R_group = MPI_COMM_WORLD;

  /* Check size of PQR_group */
  MPI_Comm_size(R_group, &size);
  if(R != size){
    if(myrank == 0){
      printf("ERROR: Incorrect number of MPI tasks. Need %d MPI tasks.", R);
    }
    error_exit(1);
  }

  if(myrank==0){
    printf("Using MPI communication layer\n");
  }
}

void
R_buf_init(int * len,
           int *  nm)
{
  /* initialize storage for nm messages in 6 directions, and associated */
  /* status information - each message is dlen long */

  double * bptr;
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
  ;

  buf_which[0] = 0;
  for(i=1; i<6; i++){
    buf_which[i] = buf_which[i-1] + which_num[i-1];
  }

  NEW( bptr, size, double * );    /* space for all messages */
  bptr_sv = bptr;   /* Save bptr for later freeing */

  NEW( buf_pool, buf_total, double ** );

  for(i = 0; i < 6; i++){    /* for each "which" set */
    k = which_len[i];       /* msg size in this set */
    for(j = 0; j < which_num[i]; j++){  /* for each msg in the set */
      buf_pool[buf_which[i] + j] = bptr;     /* ptr to next message */
      bptr = bptr + k;    /* bump over this message */
    }
    ;
  }
  ;

  NEW( buf_s_req, buf_total, MPI_Request * );
  NEW( buf_r_req, buf_total, MPI_Request * );
  NEW( buf_status, buf_total, MPI_Status * );

  buf_reset();
}

void
buf_reset()
{
  /* initialize all the things that change each iteration */
  int i;
  send_cnt = 0;

  for(i = 0; i < 6; i++){
    buf_rec[i] = 0;
  }
  for(i = 0; i < buf_total; i++){
    buf_r_req[i] = MPI_REQUEST_NULL;
    buf_s_req[i] = MPI_REQUEST_NULL;
  }
  ;     /* don't worry about status array - not used anyway */
}

void
RBufFree()
{
  /* free the memory allocated for the buffers */

  FREE( bptr_sv );
  FREE( buf_pool );
  FREE( buf_s_req );
  FREE( buf_r_req );
  FREE( buf_status );
}
