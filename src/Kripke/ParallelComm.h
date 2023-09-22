//
// Copyright (c) 2014-23, Lawrence Livermore National Security, LLC
// and Kripke project contributors. See the Kripke/COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
//

#ifndef KRIPKE_PARALLELCOMM_H__
#define KRIPKE_PARALLELCOMM_H__

#include <Kripke.h>
#include <vector>

struct Grid_Data;

namespace Kripke {

namespace Core {
  class DataStore;

  template<typename T>
  class FieldStorage;
}

class ParallelComm {
  public:
    explicit ParallelComm(Kripke::Core::DataStore &data_store);
    virtual ~ParallelComm() = default;

    // Adds a subdomain to the work queue
    virtual void addSubdomain(Kripke::Core::DataStore &data_store, SdomId sdom_id) = 0;

    // Checks if there are any outstanding subdomains to complete
    // false indicates all work is done, and all sends have completed
    virtual bool workRemaining(void);

    // Returns a vector of ready subdomains, and clears them from the ready queue
    virtual std::vector<SdomId> readySubdomains(void) = 0;

    // Marks subdomains as complete, and performs downwind communication
    virtual void markComplete(SdomId sdom_id) = 0;

  protected:
    int findSubdomain(SdomId sdom_id);
    void dequeueSubdomain(SdomId sdom_id);
    void postRecvs(Kripke::Core::DataStore &data_store, SdomId sdom_id);
    void postSends(Kripke::Core::DataStore &data_store, SdomId sdom_id_upwind, double *buffers[3]);
    void testRecieves(void);
    void waitAllSends(void);
    std::vector<SdomId> getReadyList(void);

    Kripke::Core::DataStore *m_data_store;

    Kripke::Core::FieldStorage<double> *m_plane_data[3];

    // These vectors contian the recieve requests
#ifdef KRIPKE_USE_MPI
    std::vector<MPI_Request> recv_requests;
#endif
    std::vector<int> recv_subdomains;

    // These vectors have the subdomains, and the remaining dependencies
    std::vector<int> queue_sdom_ids;
    std::vector<int> queue_depends;

    // These vectors have the remaining send requests that are incomplete
#ifdef KRIPKE_USE_MPI
    std::vector<MPI_Request> send_requests;
#endif
};


class SweepComm : public ParallelComm {
  public:
    explicit SweepComm(Kripke::Core::DataStore &data_store);
    virtual ~SweepComm();

    virtual void addSubdomain(Kripke::Core::DataStore &data_store, SdomId sdom_id);
    virtual bool workRemaining(void);
    virtual std::vector<SdomId> readySubdomains(void);
    virtual void markComplete(SdomId sdom_id);
};


class BlockJacobiComm : public ParallelComm {
  public:
    explicit BlockJacobiComm(Kripke::Core::DataStore &data_store);
    virtual ~BlockJacobiComm();

    virtual void addSubdomain(Kripke::Core::DataStore &data_store, SdomId sdom_id);
    virtual bool workRemaining(void);
    virtual std::vector<SdomId> readySubdomains(void);
    virtual void markComplete(SdomId sdom_id);

  private:
    bool posted_sends;
};


}


#endif
