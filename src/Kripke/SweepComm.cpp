#include <Kripke/SweepComm.h>
#include <Kripke/SubTVec.h>
#include <Kripke/Grid.h>

#include <fcntl.h>
#include <unistd.h>
#include <mpi.h>
#include <vector>
#include <stdio.h>


SweepComm::SweepComm(Grid_Data *data, bool bj) : grid_data(data), block_jacobi(bj)
{

}

/**
  Adds a subdomain to the work queue.
  Determines if upwind dependencies require communication, and posts appropirate Irecv's.
*/
void SweepComm::addSubdomain(int sdom_id, Subdomain &sdom){
  int mpi_rank, mpi_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  // go thru each dimensions upwind neighbors, and add the dependencies
  int num_depends = 0;
  for(int dim = 0;dim < 3;++ dim){
    // If it's a boundary condition, skip it
    if(sdom.upwind[dim].mpi_rank < 0){
      continue;
    }

    // If it's an on-rank communication (from another subdomain)
    if(sdom.upwind[dim].mpi_rank == mpi_rank){
      // skip it, but track the dependency
      num_depends ++;
      continue;
    }

    // Add request to pending list
    recv_requests.push_back(MPI_Request());
    recv_subdomains.push_back(sdom_id);

    // compute the tag id of THIS subdomain (tags are always based on destination)
    int tag = mpi_rank + mpi_size*sdom_id;

    // Post the recieve
    MPI_Irecv(sdom.plane_data[dim]->ptr(), sdom.plane_data[dim]->elements, MPI_DOUBLE, sdom.upwind[dim].mpi_rank,
      tag, MPI_COMM_WORLD, &recv_requests[recv_requests.size()-1]);

    // increment number of dependencies
    num_depends ++;
  }

  // add subdomain to queue
  queue_sdom_ids.push_back(sdom_id);
  queue_subdomains.push_back(&sdom);
  queue_depends.push_back(num_depends);


}

// Checks if there are any outstanding subdomains to complete
// false indicates all work is done, and all sends have completed
bool SweepComm::workRemaining(void){
  // If there are outstanding subdomains to process, return true
  if(recv_requests.size() > 0 || queue_subdomains.size() > 0){
    return true;
  }

  // Wait for all remaining sends to complete, then return false
  int num_sends = send_requests.size();
  if(num_sends > 0){
    std::vector<MPI_Status> status(num_sends);
    MPI_Waitall(num_sends, &send_requests[0], &status[0]);
    send_requests.clear();
  }

  return false;
}

/**
  Checks for incomming messages, and returns a list of ready subdomain id's
*/
std::vector<int> SweepComm::readySubdomains(void){
  // If we're in BJ mode, theo first call to readySubdomains will perform
  // all sends, wait for sends to complete, then report everything as
  // ready
  if(block_jacobi){
    // have we not sent everything?
    if(recv_requests.size() > 0){
      // loop over subdomains and send
      for(int i = 0;i < queue_subdomains.size();++ i){
        postSends(queue_subdomains[i]);
      }

      // wait for all recieves to complete
      int num_recvs = recv_requests.size();
      std::vector<MPI_Status> status(num_recvs);
      MPI_Waitall(num_recvs, &recv_requests[0], &status[0]);
      recv_requests.clear();
      recv_subdomains.clear();
    }


    // return list of all subdomains in queue
    return queue_sdom_ids;
  }


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
      for(int i = 0;i < queue_sdom_ids.size();++ i){
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


  // build up a list of ready subdomains
  std::vector<int> ready;
  for(int i = 0;i < queue_depends.size();++ i){
    if(queue_depends[i] == 0){
      ready.push_back(queue_sdom_ids[i]);
    }
  }
  return ready;
}

/**
  Marks subdomains as complete, and performs downwind communication
*/
int SweepComm::findSubdomain(int sdom_id){

  // find subdomain in queue
  int index;
  for(index = 0;index < queue_sdom_ids.size();++ index){
    if(queue_sdom_ids[index] == sdom_id){
      break;
    }
  }
  if(index == queue_sdom_ids.size()){
    printf("Cannot find subdomain id %d in work queue\n", sdom_id);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  return index;
}

void SweepComm::markComplete(int sdom_id){
  int index = findSubdomain(sdom_id);

  // Get subdomain pointer before removing it from queue
  Subdomain *sdom = queue_subdomains[index];

  // remove subdomain from queue
  queue_sdom_ids.erase(queue_sdom_ids.begin()+index);
  queue_subdomains.erase(queue_subdomains.begin()+index);
  queue_depends.erase(queue_depends.begin()+index);

  // Send downwind info if we are doing a normal sweep
  if(!block_jacobi){
    postSends(sdom);
  }
}

void SweepComm::postSends(Subdomain *sdom){
  // post sends for downwind dependencies
  int mpi_rank, mpi_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  for(int dim = 0;dim < 3;++ dim){
    // If it's a boundary condition, skip it
    if(sdom->downwind[dim].mpi_rank < 0){
      continue;
    }

    // If it's an on-rank communication (to another subdomain)
    if(sdom->downwind[dim].mpi_rank == mpi_rank){
      // find the local subdomain in the queue, and decrement the counter
      for(int i = 0;i < queue_sdom_ids.size();++ i){
        if(queue_sdom_ids[i] == sdom->downwind[dim].subdomain_id){
          queue_depends[i] --;
          break;
        }
      }

      // copy the boundary condition data into the downwinds plane data
      Subdomain &sdom_downwind = grid_data->subdomains[sdom->downwind[dim].subdomain_id];
      sdom_downwind.plane_data[dim]->copy(*sdom->plane_data[dim]);
      int num_elem = sdom_downwind.plane_data[dim]->elements;
      double const * KRESTRICT src_ptr = sdom->plane_data[dim]->ptr();
      double * KRESTRICT dst_ptr = sdom_downwind.plane_data[dim]->ptr();
      for(int i = 0;i < num_elem;++ i){
        dst_ptr[i] = src_ptr[i];
      }

      continue;
    }

    // At this point, we know that we have to send an MPI message
    // Add request to send queue
    send_requests.push_back(MPI_Request());

    // compute the tag id of TARGET subdomain (tags are always based on destination)
    int tag = sdom->downwind[dim].mpi_rank + mpi_size*sdom->downwind[dim].subdomain_id;

    // Post the send
    MPI_Isend(sdom->plane_data[dim]->ptr(), sdom->plane_data[dim]->elements, MPI_DOUBLE, sdom->downwind[dim].mpi_rank,
      tag, MPI_COMM_WORLD, &send_requests[send_requests.size()-1]);
  }
}

