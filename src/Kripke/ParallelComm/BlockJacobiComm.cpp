/*
 * NOTICE
 *
 * This work was produced at the Lawrence Livermore National Laboratory (LLNL)
 * under contract no. DE-AC-52-07NA27344 (Contract 44) between the U.S.
 * Department of Energy (DOE) and Lawrence Livermore National Security, LLC
 * (LLNS) for the operation of LLNL. The rights of the Federal Government are
 * reserved under Contract 44.
 *
 * DISCLAIMER
 *
 * This work was prepared as an account of work sponsored by an agency of the
 * United States Government. Neither the United States Government nor Lawrence
 * Livermore National Security, LLC nor any of their employees, makes any
 * warranty, express or implied, or assumes any liability or responsibility
 * for the accuracy, completeness, or usefulness of any information, apparatus,
 * product, or process disclosed, or represents that its use would not infringe
 * privately-owned rights. Reference herein to any specific commercial products,
 * process, or service by trade name, trademark, manufacturer or otherwise does
 * not necessarily constitute or imply its endorsement, recommendation, or
 * favoring by the United States Government or Lawrence Livermore National
 * Security, LLC. The views and opinions of authors expressed herein do not
 * necessarily state or reflect those of the United States Government or
 * Lawrence Livermore National Security, LLC, and shall not be used for
 * advertising or product endorsement purposes.
 *
 * NOTIFICATION OF COMMERCIAL USE
 *
 * Commercialization of this product is prohibited without notifying the
 * Department of Energy (DOE) or Lawrence Livermore National Security.
 */

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
  Set const &set_iplane = data_store.getVariable<Set>("Set/IPlane");
  Set const &set_jplane = data_store.getVariable<Set>("Set/JPlane");
  Set const &set_kplane = data_store.getVariable<Set>("Set/KPlane");
  data_store.newVariable<Field_IPlane>("old_i_plane", set_iplane);
  data_store.newVariable<Field_JPlane>("old_j_plane", set_jplane);
  data_store.newVariable<Field_KPlane>("old_k_plane", set_kplane);
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


