//
// Copyright (c) 2014-19, Lawrence Livermore National Security, LLC
// and Kripke project contributors. See the COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
//

#include <Kripke/SweepSolver.h>

#include <Kripke.h>
#include <Kripke/Kernel.h>
#include <Kripke/ParallelComm.h>
#include <Kripke/Timing.h>
#include <Kripke/VarTypes.h>
#include <vector>
#include <stdio.h>

using namespace Kripke;

/**
  Perform full parallel sweep algorithm on subset of subdomains.
*/
void Kripke::SweepSolver (Kripke::Core::DataStore &data_store, std::vector<SdomId> subdomain_list, bool block_jacobi)
{
  KRIPKE_TIMER(data_store, SweepSolver);

  // Initialize plane data
  Kripke::Kernel::kConst(data_store.getVariable<Field_IPlane>("i_plane"), 0.0);
  Kripke::Kernel::kConst(data_store.getVariable<Field_JPlane>("j_plane"), 0.0);
  Kripke::Kernel::kConst(data_store.getVariable<Field_KPlane>("k_plane"), 0.0);

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
//    Kripke::Core::Comm default_comm;
//    printf("SweepSolver: rank=%d, sdom=%d\n", (int)default_comm.rank(), (int)*subdomain_list[i]);
    SdomId sdom_id = subdomain_list[i];
    comm->addSubdomain(data_store, sdom_id);
  }

  auto &field_upwind = data_store.getVariable<Field_Adjacency>("upwind");

  /* Loop until we have finished all of our work */
  while(comm->workRemaining()){

    std::vector<SdomId> sdom_ready = comm->readySubdomains();
    int backlog = sdom_ready.size();

    // Run top of list
    if(backlog > 0){
      SdomId sdom_id = sdom_ready[0];

      auto upwind = field_upwind.getView(sdom_id);

      // Clear boundary conditions
      if(upwind(Dimension{0}) == -1){
        Kripke::Kernel::kConst(data_store.getVariable<Field_IPlane>("i_plane"), sdom_id, 0.0);
      }
      if(upwind(Dimension{1}) == -1){
        Kripke::Kernel::kConst(data_store.getVariable<Field_JPlane>("j_plane"), sdom_id, 0.0);
      }
      if(upwind(Dimension{2}) == -1){
        Kripke::Kernel::kConst(data_store.getVariable<Field_KPlane>("k_plane"), sdom_id, 0.0);
      }

      // Perform subdomain sweep
      Kripke::Kernel::sweepSubdomain(data_store, Kripke::SdomId{sdom_id});

      // Mark as complete (and do any communication)
      comm->markComplete(sdom_id);
    }
  }

  delete comm;

//  printf("\nAfter sweep psi:\n");
//  data_store.getVariable<Field_Flux>("psi").dump();

}


