/* Common declarations for the functions in comm.c */
#ifndef KRIPKE_COMM_H__
#define KRIPKE_COMM_H__

#include<vector>
#include<mpi.h>

class Grid_Data;
class Subdomain;

class SweepComm {
  public:
    SweepComm(Grid_Data *data, bool bj_mode = false);

    // Adds a subdomain to the work queue
    void addSubdomain(int sdom_id, Subdomain &sdom);

    // Checks if there are any outstanding subdomains to complete
    // false indicates all work is done, and all sends have completed
    bool workRemaining(void);

    // Returns a vector of ready subdomains, and clears them from the ready queue
    std::vector<int> readySubdomains(void);

    // Marks subdomains as complete, and performs downwind communication
    void markComplete(int sdom_id);

  private:
    int findSubdomain(int sdom_id);
    void postSends(Subdomain *sdom, int index);

    Grid_Data *grid_data;
    bool block_jacobi; // Turns off parallel sweep, does up-front boundary exchange

    // These vectors contian the recieve requests
    std::vector<MPI_Request> recv_requests;
    std::vector<int> recv_subdomains;

    // These vectors have the subdomains, and the remaining dependencies
    std::vector<int> queue_sdom_ids;
    std::vector<Subdomain *> queue_subdomains;
    std::vector<int> queue_depends;
    std::vector<std::vector<double> > bj_send_buffers[3];

    // These vectors have the remaining send requests that are incomplete
    std::vector<MPI_Request> send_requests;
};



#endif
