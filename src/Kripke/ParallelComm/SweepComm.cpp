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

