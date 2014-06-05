/* Common declarations for the functions in comm.c */
#ifndef KRIPKE_COMM_H__
#define KRIPKE_COMM_H__

#include<vector>
#include<mpi.h>


/**
 * This class provides the asynchronous point-to-point message passing needed
 * for the MPI-based sweep algorithm in Sweep_Solver.cpp
 *
 */
class Comm {
  public:
    Comm(int * len, int * nm);

    void R_recv_dir (int which, int r_member );
    int  R_recv_test (int which, double **msg );
    void R_send (double * msg, int r_member, int length );
    void R_wait_send();

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


#endif
