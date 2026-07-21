//
// Copyright (c) 2014-25, Lawrence Livermore National Security, LLC
// and Kripke project contributors. See the Kripke/COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
//

#include <Kripke/ParallelComm.h>

#include <Kripke/Core/Comm.h>
#include <Kripke/Core/Field.h>
#include <Kripke/VarTypes.h>

using namespace Kripke;

namespace {

bool useGpuAwareMPI(Kripke::Core::FieldStorage<double> &field)
{
#ifdef KRIPKE_USE_GPU_AWARE_MPI
  return field.getAllocationSpace() == chai::GPU;
#else
  (void)field;
  return false;
#endif
}

void synchronizeDeviceForMPI(void)
{
#if defined(KRIPKE_USE_CUDA)
  RAJA::synchronize<RAJA::cuda_synchronize>();
#elif defined(KRIPKE_USE_HIP)
  RAJA::synchronize<RAJA::hip_synchronize>();
#endif
}

} // namespace

// Helper for copying plane data between two GPU allocations, only used in ParallelComm
static void copyPlane(Kripke::Core::FieldStorage<double> &dst_plane,
                      Kripke::SdomId dst_sdom_id,
                      Kripke::Core::FieldStorage<double> &src_plane,
                      Kripke::SdomId src_sdom_id)
{
  int num_elem = src_plane.size(src_sdom_id);
  KRIPKE_ASSERT(dst_plane.size(dst_sdom_id) == (size_t)num_elem,
      "Cannot copy plane data with different subdomain sizes");

#if defined(KRIPKE_USE_CHAI) && (defined(KRIPKE_USE_CUDA) || defined(KRIPKE_USE_HIP))
  if(dst_plane.getAllocationSpace() == chai::GPU &&
     src_plane.getAllocationSpace() == chai::GPU){
    double *dst = dst_plane.getDeviceData(dst_sdom_id);
    double const *src = src_plane.getDeviceData(src_sdom_id);
#if defined(KRIPKE_USE_CUDA)
    using PlaneCopyExec = RAJA::cuda_exec<256>;
#else
    using PlaneCopyExec = RAJA::hip_exec<256>;
#endif
    RAJA::forall<PlaneCopyExec>(
      RAJA::RangeSegment(0, num_elem),
      KRIPKE_LAMBDA (RAJA::Index_type i){
        dst[i] = src[i];
    });
    return;
  }
#endif  // #if defined(KRIPKE_USE_CHAI) && (defined(KRIPKE_USE_CUDA) || defined(KRIPKE_USE_HIP))

  // Fallback for non-CHAI and CHAI fields that are not both GPU-backed.
  double *dst = dst_plane.getHostData(dst_sdom_id);
  double const *src = src_plane.getHostDataConst(src_sdom_id);
  for(int i = 0;i < num_elem;++ i){
    dst[i] = src[i];
  }
}

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
  Receives use either direct GPU plane buffers or direct host plane buffers.
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
    recv_dimensions.push_back(*dim);

    auto &plane_data = *m_plane_data[*dim];
    size_t plane_data_size = plane_data.size(sdom_id);
    double *plane_data_ptr = nullptr;

    if(useGpuAwareMPI(plane_data)){
      plane_data_ptr = plane_data.getDeviceData(sdom_id);
    }
    else{
      plane_data_ptr = plane_data.getHostData(sdom_id);
    }

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
                             Kripke::Core::FieldStorage<double> *src_plane_data[3])
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

      // copy the boundary condition data into the downwind plane data
      copyPlane(*m_plane_data[*dim], sdom_id_downwind,
                *src_plane_data[*dim], sdom_id);
      continue;
    }

#ifdef KRIPKE_USE_MPI

    // At this point, we know that we have to send an MPI message
    // Add request to send queue
    send_requests.push_back(MPI_Request());

    // Get size of outgoing boudnary data
    auto &src_plane = *src_plane_data[*dim];
    size_t plane_data_size = src_plane.size(sdom_id);
    double *src_buffer = nullptr;

    if(useGpuAwareMPI(src_plane)){
      synchronizeDeviceForMPI();
      src_buffer = src_plane.getDeviceData(sdom_id);
    }
    else{
      src_buffer = const_cast<double *>(src_plane.getHostDataConst(sdom_id));
    }

    // Post the send
    MPI_Isend(src_buffer, plane_data_size, MPI_DOUBLE, downwind_rank,
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

#ifdef KRIPKE_USE_GPU_AWARE_MPI
      if(useGpuAwareMPI(*m_plane_data[recv_dimensions[index]])){
        m_plane_data[recv_dimensions[index]]->registerDeviceTouch(SdomId{sdom_id});
      }
#endif

      // remove the request from the list
      recv_requests.erase(recv_requests.begin()+index);
      recv_subdomains.erase(recv_subdomains.begin()+index);
      recv_dimensions.erase(recv_dimensions.begin()+index);
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
