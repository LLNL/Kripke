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

#include <Kripke/SweepSolver.h>
#include <Kripke.h>
#include <Kripke/Grid.h>
#include <Kripke/ParallelComm.h>
#include <Kripke/Subdomain.h>
#include <Kripke/SubTVec.h>
#include <Kripke/VarTypes.h>
#include <vector>
#include <stdio.h>

using namespace Kripke;

/**
  Perform full parallel sweep algorithm on subset of subdomains.
*/
void Kripke::SweepSolver (Kripke::DataStore &data_store, std::vector<SdomId> subdomain_list, Grid_Data *grid_data, bool block_jacobi)
{
  KRIPKE_TIMER(data_store, SweepSolver);

  // Create a new sweep communicator object
  ParallelComm *comm = NULL;
  if(block_jacobi){
    comm = new BlockJacobiComm(data_store);
  }
  else {
    comm = new SweepComm(data_store);
  }

  // Add all subdomains in our list
  for(size_t i = 0;i < subdomain_list.size();++ i){
    SdomId sdom_id = subdomain_list[i];
    comm->addSubdomain(sdom_id, grid_data->subdomains[*sdom_id]);
  }


  /* Loop until we have finished all of our work */
  while(comm->workRemaining()){

    std::vector<SdomId> sdom_ready = comm->readySubdomains();
    int backlog = sdom_ready.size();

    // Run top of list
    if(backlog > 0){
      SdomId sdom_id = sdom_ready[0];
      Subdomain &sdom = grid_data->subdomains[*sdom_id];
      // Clear boundary conditions
      if(sdom.upwind[0].subdomain_id == -1){
        Kripke::Kernel::kConst(data_store.getVariable<Field_IPlane>("i_plane"), 0.0);
      }
      if(sdom.upwind[1].subdomain_id == -1){
        Kripke::Kernel::kConst(data_store.getVariable<Field_JPlane>("j_plane"), 0.0);
      }
      if(sdom.upwind[2].subdomain_id == -1){
        Kripke::Kernel::kConst(data_store.getVariable<Field_KPlane>("k_plane"), 0.0);
      }

      // Perform subdomain sweep
      Kripke::Kernel::sweepSubdomain(data_store, Kripke::SdomId{sdom_id});

      // Mark as complete (and do any communication)
      comm->markComplete(sdom_id);
    }
  }

  delete comm;
}


