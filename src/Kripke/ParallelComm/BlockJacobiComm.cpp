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
using namespace Kripke::Core;

BlockJacobiComm::BlockJacobiComm(Kripke::Core::DataStore &data_store) :
ParallelComm(data_store), posted_sends(false)
{
  ArchLayoutV al_v = data_store.getVariable<ArchLayout>("al").al_v;

  Set const &set_iplane = data_store.getVariable<Set>("Set/IPlane");
  Set const &set_jplane = data_store.getVariable<Set>("Set/JPlane");
  Set const &set_kplane = data_store.getVariable<Set>("Set/KPlane");
  createField<Field_IPlane>(data_store, "old_i_plane", al_v, set_iplane);
  createField<Field_JPlane>(data_store, "old_j_plane", al_v, set_jplane);
  createField<Field_KPlane>(data_store, "old_k_plane", al_v, set_kplane);
}

BlockJacobiComm::~BlockJacobiComm(){
  m_data_store->deleteVariable("old_i_plane");
  m_data_store->deleteVariable("old_j_plane");
  m_data_store->deleteVariable("old_k_plane");
}

/**
  Adds a subdomain to the work queue.
  Determines if upwind dependencies require communication, and posts appropirate Irecv's.
*/
void BlockJacobiComm::addSubdomain(Kripke::Core::DataStore &data_store, SdomId sdom_id){

  // Copy old flux data to send buffers
  auto &i_plane = m_data_store->getVariable<Field_IPlane>("i_plane");
  auto &old_i_plane = m_data_store->getVariable<Field_IPlane>("old_i_plane");
  Kernel::kCopy(old_i_plane, i_plane);

  auto &j_plane = m_data_store->getVariable<Field_JPlane>("j_plane");
  auto &old_j_plane = m_data_store->getVariable<Field_JPlane>("old_j_plane");
  Kernel::kCopy(old_j_plane, j_plane);

  auto &k_plane = m_data_store->getVariable<Field_KPlane>("k_plane");
  auto &old_k_plane = m_data_store->getVariable<Field_KPlane>("old_k_plane");
  Kernel::kCopy(old_k_plane, k_plane);

  // post recieves
  postRecvs(data_store, sdom_id);

}

// Checks if there are any outstanding subdomains to complete
// false indicates all work is done, and all sends have completed
bool BlockJacobiComm::workRemaining(void){
  if(!posted_sends){

    auto &old_i_plane = m_data_store->getVariable<Field_IPlane>("old_i_plane");
    auto &old_j_plane = m_data_store->getVariable<Field_JPlane>("old_j_plane");
    auto &old_k_plane = m_data_store->getVariable<Field_KPlane>("old_k_plane");

    // post sends for all queued subdomains
    for(size_t i = 0;i < queue_sdom_ids.size();++ i){
      SdomId sdom_id(queue_sdom_ids[i]);

      // Send new downwind info for sweep
      double *buf[3] = {
          old_i_plane.getData(sdom_id),
          old_j_plane.getData(sdom_id),
          old_k_plane.getData(sdom_id)
      };

      postSends(*m_data_store, sdom_id, buf);
    }
    posted_sends = true;
  }
  // Since we communicate fluxes before local sweeps, when we are
  // out of work, there is no further synchronization
  if(ParallelComm::workRemaining()){
    return true;
  }
  waitAllSends();

  return false;
}

/**
  Checks for incomming messages, and returns a list of ready subdomain id's
*/
std::vector<SdomId> BlockJacobiComm::readySubdomains(void){
  testRecieves();

  // return list of any ready subdomains
  return getReadyList();
}



void BlockJacobiComm::markComplete(SdomId sdom_id){
  // remove subdomain from work queue
  dequeueSubdomain(sdom_id);
}


