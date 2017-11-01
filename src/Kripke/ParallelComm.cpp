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

#include <Kripke/Core/Comm.h>
#include <Kripke/Core/Field.h>
#include <Kripke/VarTypes.h>

using namespace Kripke;

ParallelComm::ParallelComm(Kripke::Core::DataStore &data_store) :
  m_data_store(&data_store)
{
  m_plane_data[0] = &m_data_store->getVariable<Field_IPlane>("i_plane");
  m_plane_data[1] = &m_data_store->getVariable<Field_JPlane>("j_plane");
  m_plane_data[2] = &m_data_store->getVariable<Field_KPlane>("k_plane");

}



/**
  Finds subdomain in the queue by its subdomain id.
*/
int ParallelComm::findSubdomain(SdomId sdom_id){

  // find subdomain in queue
  size_t index;
  for(index = 0;index < queue_sdom_ids.size();++ index){
    if(queue_sdom_ids[index] == *sdom_id){
      break;
    }
  }
  if(index == queue_sdom_ids.size()){
    KRIPKE_ABORT("Cannot find subdomain id %ld in work queue\n", (long)*sdom_id);
  }

  return index;
}


void ParallelComm::dequeueSubdomain(SdomId sdom_id){
  int index = findSubdomain(sdom_id);

  // remove subdomain from queue
  queue_sdom_ids.erase(queue_sdom_ids.begin()+index);
  queue_depends.erase(queue_depends.begin()+index);

}

/**
  Adds a subdomain to the work queue.
  Determines if upwind dependencies require communication, and posts appropriate Irecv's.
  All receives use the plane_data[] arrays as receive buffers.
*/
void ParallelComm::postRecvs(Kripke::Core::DataStore &data_store, SdomId sdom_id){
  using namespace Kripke::Core;
  Comm comm;
  int mpi_rank = comm.rank();

  auto upwind = data_store.getVariable<Field_Adjacency>("upwind").getView(sdom_id);

  auto global_to_rank = data_store.getVariable<Field_GlobalSdomId2Rank>("GlobalSdomId2Rank").getView(SdomId{0});

#ifdef KRIPKE_USE_MPI
  auto local_to_global = data_store.getVariable<Field_SdomId2GlobalSdomId>("SdomId2GlobalSdomId").getView(SdomId{0});
#endif

  // go thru each dimensions upwind neighbors, and add the dependencies
  int num_depends = 0;
  for(Dimension dim{0};dim < 3;++ dim){

    // If it's a boundary condition, skip it
    if(upwind(dim) < 0){
      continue;
    }

    // If it's an on-rank communication (from another subdomain)
    GlobalSdomId upwind_sdom = upwind(dim);
    int upwind_rank = global_to_rank(upwind_sdom);

    if(upwind_rank == mpi_rank){
      // skip it, but track the dependency
      num_depends ++;
      continue;
    }

#ifdef KRIPKE_USE_MPI

    // Add request to pending list
    recv_requests.push_back(MPI_Request());
    recv_subdomains.push_back(*sdom_id);

    auto &plane_data = *m_plane_data[*dim];
    double *plane_data_ptr = plane_data.getData(sdom_id);
    size_t plane_data_size = plane_data.size(sdom_id);

    GlobalSdomId global_sdom_id = local_to_global(sdom_id);

    // Post the recieve
    MPI_Irecv(plane_data_ptr, plane_data_size, MPI_DOUBLE, upwind_rank,
      *global_sdom_id, MPI_COMM_WORLD, &recv_requests[recv_requests.size()-1]);

    // increment number of dependencies
    num_depends ++;
#else
    // No MPI, so this doesn't make sense
    KRIPKE_ASSERT("All comms should be on-node without MPI");
#endif
  }

  // add subdomain to queue
  queue_sdom_ids.push_back(*sdom_id);
  queue_depends.push_back(num_depends);
}

void ParallelComm::postSends(Kripke::Core::DataStore &data_store, Kripke::SdomId sdom_id,
                             double *src_buffers[3])
{
  // post sends for downwind dependencies
  Kripke::Core::Comm comm;
  int mpi_rank = comm.rank();

  auto downwind = data_store.getVariable<Field_Adjacency>("downwind").getView(sdom_id);
  auto global_to_rank = data_store.getVariable<Field_GlobalSdomId2Rank>("GlobalSdomId2Rank").getView(SdomId{0});
  auto global_to_sdom_id = data_store.getVariable<Field_GlobalSdomId2SdomId>("GlobalSdomId2SdomId").getView(SdomId{0});

  for(Dimension dim{0};dim < 3;++ dim){
    // If it's a boundary condition, skip it
    if(downwind(dim) < 0){
      continue;
    }



    // If it's an on-rank communication (to another subdomain)
    GlobalSdomId downwind_sdom = downwind(dim);
    int downwind_rank = global_to_rank(downwind_sdom);
    if(downwind_rank == mpi_rank){

      SdomId sdom_id_downwind = global_to_sdom_id(downwind(dim));

      // find the local subdomain in the queue, and decrement the counter
      for(size_t i = 0;i < queue_sdom_ids.size();++ i){
        if(queue_sdom_ids[i] == *sdom_id_downwind){
          queue_depends[i] --;
          break;
        }
      }

      // copy the boundary condition data into the downwinds plane data
      auto dst_plane = m_plane_data[*dim]->getView1d(sdom_id_downwind);
      int num_elem = m_plane_data[*dim]->size(sdom_id_downwind);
      for(int i = 0;i < num_elem;++ i){
        dst_plane(i) = src_buffers[*dim][i];
      }
      continue;
    }

#ifdef KRIPKE_USE_MPI

    // At this point, we know that we have to send an MPI message
    // Add request to send queue
    send_requests.push_back(MPI_Request());

    // Get size of outgoing boudnary data
    auto &plane_data = *m_plane_data[*dim];
    size_t plane_data_size = plane_data.size(sdom_id);

    // Post the send
    MPI_Isend(src_buffers[*dim], plane_data_size, MPI_DOUBLE, downwind_rank,
      *downwind_sdom, MPI_COMM_WORLD, &send_requests[send_requests.size()-1]);

#else
    // We cannot SEND anything without MPI, so fail
    KRIPKE_ASSERT("Cannot send messages without MPI");
#endif

  }
}


// Checks if there are any outstanding subdomains to complete
bool ParallelComm::workRemaining(void){
#ifdef KRIPKE_USE_MPI
  return (recv_requests.size() > 0 || queue_sdom_ids.size() > 0);
#else
  return (queue_sdom_ids.size() > 0);
#endif
}


// Blocks until all sends have completed, and flushes the send queues
void ParallelComm::waitAllSends(void){
#ifdef KRIPKE_USE_MPI
  // Wait for all remaining sends to complete, then return false
  int num_sends = send_requests.size();
  if(num_sends > 0){
    std::vector<MPI_Status> status(num_sends);
    MPI_Waitall(num_sends, &send_requests[0], &status[0]);
    send_requests.clear();
  }
#endif
}

/**
  Checks for incomming messages, and does relevant bookkeeping.
*/
void ParallelComm::testRecieves(void){
#ifdef KRIPKE_USE_MPI
  // Check for any recv requests that have completed
  int num_requests = recv_requests.size();
  bool done = false;
  while(!done && num_requests > 0){
    // Create array of status variables
    std::vector<MPI_Status> recv_status(num_requests);

    // Ask if either one or none of the recvs have completed?
    int index; // this will be the index of request that completed
    int complete_flag; // this is set to TRUE if somthing completed
    MPI_Testany(num_requests, &recv_requests[0], &index, &complete_flag, &recv_status[0]);

    if(complete_flag != 0){

      // get subdomain that this completed for
      int sdom_id = recv_subdomains[index];

      // remove the request from the list
      recv_requests.erase(recv_requests.begin()+index);
      recv_subdomains.erase(recv_subdomains.begin()+index);
      num_requests --;

      // decrement the dependency count for that subdomain
      for(size_t i = 0;i < queue_sdom_ids.size();++ i){
        if(queue_sdom_ids[i] == sdom_id){
          queue_depends[i] --;
          break;
        }
      }
    }
    else{
      done = true;
    }
  }
#endif
}


std::vector<SdomId> ParallelComm::getReadyList(void){
  // build up a list of ready subdomains
  std::vector<SdomId> ready;
  for(size_t i = 0;i < queue_depends.size();++ i){
    if(queue_depends[i] == 0){
      ready.push_back(SdomId(queue_sdom_ids[i]));
    }
  }
  return ready;
}
