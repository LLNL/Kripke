/* Common declarations for the functions in comm.c */
#ifndef KRIPKE_COMM_H__
#define KRIPKE_COMM_H__

#include<mpi.h>

/*================= POINT TO POINT MESSAGE PASSING ==================*/

extern void R_buf_init (int * len, int * nm);
/*
   R_buf_init() initializes storage for messages (and associated status arrays)
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

extern void RBufFree ();
/*
   frees all the memory allocated for the buffers used in sweeping.
*/

extern void R_recv_dir (int which, int r_member );
/*
   R_recv_dir() posts receives in the signed direction ``which''.  The number
   of receives posted is nm[which], where nm[] is the array passed to
   the initialization function R_buf_init().  The argument ``r_member''
   is the r coordinate of the node from which the nm[which] messages will
   be sent.  R_recv_dir is non-blocking, i.e., it returns after posting
   the receives without waiting for any messages to arrive.
*/

extern int  R_recv_test (int which, double **msg );
/*
   R_recv_test() checks to see if any of the posted receives in the
   signed direction ``which'' has been satisfied (i.e, a message has
   arrived).  If no received message is pending, the function returns 0.
   If a message has arrived, its address is assigned to *msg and the
   function returns 1.
*/

extern void R_send (double * msg, int r_member, int length );
/*
   Called from a node (p,q,r), R_send() sends the message of ``length''
   doubles pointed at by ``msg'' to node (p,q,r_member).  The send is
   non-blocking, i.e., it returns immediately after submitting the message
   to the underlying communication system.
*/

extern void R_wait_send();
/*
   R_wait_send() waits for all posted sends to complete.
*/

/*================= Miscellaneous Functions ==================*/

extern void error_exit( int flag );
/*
   error_exit() aborts the concurrent execution, passing the flag value
   to lower-level (system-dependent) abort routines, which may ignore it.
*/

extern MPI_Comm GetRGroup( );
/*
   Returns the R_group communicator.
*/

extern int GetRrank( );
/*
   Returns the R_group rank.
*/

extern void create_R_grid( int R );
/*
   create_PQR() creates the R grid.
*/

#endif
