//
// Copyright (c) 2014-19, Lawrence Livermore National Security, LLC
// and Kripke project contributors. See the COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
//

#include <Kripke/ParallelComm.h>

#include <Kripke/Core/Field.h>
#include <Kripke/Kernel.h>
#include <Kripke/Timing.h>
#include <Kripke/VarTypes.h>

#include <fcntl.h>
#include <unistd.h>
#include <vector>
#include <stdio.h>

using namespace Kripke;


SweepComm::SweepComm(Kripke::Core::DataStore &data_store) : ParallelComm(data_store)
{

}

SweepComm::~SweepComm(){
}

/**
  Adds a subdomain to the work queue.
  Determines if upwind dependencies require communication, and posts appropirate Irecv's.
*/
void SweepComm::addSubdomain(Kripke::Core::DataStore &data_store, SdomId sdom_id){
  // Post recieves for upwind dependencies, and add to the queue
  postRecvs(data_store, sdom_id);
}


// Checks if there are any outstanding subdomains to complete
// false indicates all work is done, and all sends have completed
bool SweepComm::workRemaining(void){
  // If there are outstanding subdomains to process, return true
  if(ParallelComm::workRemaining()){
    return true;
  }

  // No more work, so make sure all of our sends have completed
  // before we continue
  waitAllSends();

  return false;
}


/**
  Checks for incomming messages, and returns a list of ready subdomain id's
*/
std::vector<SdomId> SweepComm::readySubdomains(void){
  // check for incomming messages
  testRecieves();

  // build up a list of ready subdomains
  return getReadyList();
}


void SweepComm::markComplete(SdomId sdom_id){
  // Get subdomain pointer and remove from work queue
  dequeueSubdomain(sdom_id);

  auto &i_plane = m_data_store->getVariable<Field_IPlane>("i_plane");
  auto &j_plane = m_data_store->getVariable<Field_JPlane>("j_plane");
  auto &k_plane = m_data_store->getVariable<Field_KPlane>("k_plane");

  // Send new downwind info for sweep
  double *buf[3] = {
    i_plane.getData(sdom_id),
    j_plane.getData(sdom_id),
    k_plane.getData(sdom_id)
  };
  postSends(*m_data_store, sdom_id, buf);
}

