/* Common declarations for the functions in comm.c */
#ifndef KRIPKE_COMM_H__
#define KRIPKE_COMM_H__

#include<vector>
#include<mpi.h>

/*================= POINT TO POINT MESSAGE PASSING ==================*/

class Comm {
  public:
    Comm(int * len, int * nm);
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

    ~Comm();
    /*
       frees all the memory allocated for the buffers used in sweeping.
    */

    void R_recv_dir (int which, int r_member );
    /*
       R_recv_dir() posts receives in the signed direction ``which''.  The number
       of receives posted is nm[which], where nm[] is the array passed to
       the initialization function R_buf_init().  The argument ``r_member''
       is the r coordinate of the node from which the nm[which] messages will
       be sent.  R_recv_dir is non-blocking, i.e., it returns after posting
       the receives without waiting for any messages to arrive.
    */

    int  R_recv_test (int which, double **msg );
    /*
       R_recv_test() checks to see if any of the posted receives in the
       signed direction ``which'' has been satisfied (i.e, a message has
       arrived).  If no received message is pending, the function returns 0.
       If a message has arrived, its address is assigned to *msg and the
       function returns 1.
    */

    void R_send (double * msg, int r_member, int length );
    /*
       Called from a node (p,q,r), R_send() sends the message of ``length''
       doubles pointed at by ``msg'' to node (p,q,r_member).  The send is
       non-blocking, i.e., it returns immediately after submitting the message
       to the underlying communication system.
    */

    void R_wait_send();
    /*
       R_wait_send() waits for all posted sends to complete.
    */

    /*================= Miscellaneous Functions ==================*/




    MPI_Comm GetRGroup( );
    /*
       Returns the R_group communicator.
    */

    int GetRrank( );
    /*
       Returns the R_group rank.
    */

    void create_R_grid( int R );
    /*
       create_PQR() creates the R grid.
    */

  private:
    void buf_reset();

    /*  -------------------------- message buffer declarations -------------- */
    /* message buffer information, used in send-receive calls */
    /* the length of several arrays is not known until runtime, */
    /* so will be allocated via malloc in the buf_init routine   */

    /* these items are initialized in buf_init and don't change */
    int buf_total;     /* count of total number of messages */
    int which_num[6];  /* max number of message in each direction */
    int which_len[6];  /* size of messages in each direction */
    int buf_which[6];  /* array to say where in pool each direction starts */
    std::vector<double *> buf_pool; /* array of  pointers to each buffer in pool */
    std::vector<double> bptr_sv;  /* pointer to memory for all buffers */

    /* these items are reinitialized every iteration through the main code */
    int send_cnt;     /* index into buf_s_req for next send-handle */
    int buf_rec[6];   /* array to track whether which-direction is boundary */
                      /* if value is 1 - receives have been posted in direction */
                      /* if values is > 1 error: rcv's posted more than once*/
                      /* if value is <= 0 direction may be used for boundry bufs*/
                      /* amount less than  0 tracks which to allocate next  */
    std::vector<MPI_Request> buf_r_req; /* array to hold async recv message handles */
                            /* entries 1-1 with buf_pool, initialed to null-handle
                              */
    std::vector<MPI_Request> buf_s_req; /* array to hold async send message handles */
                            /* entries made in order send's posted (see send_cnt)
                              */
    std::vector<MPI_Status> buf_status; /* not used except for debug and MPI returns it */

};

void error_exit( int flag );
/*
     error_exit() aborts the concurrent execution, passing the flag value
     to lower-level (system-dependent) abort routines, which may ignore it.
  */

#endif
