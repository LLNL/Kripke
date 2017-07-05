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

#include <Kripke/SteadyStateSolver.h>
#include <Kripke.h>
#include <Kripke/Comm.h>
#include <Kripke/Grid.h>
#include <Kripke/ParallelComm.h>
#include <Kripke/Subdomain.h>
#include <Kripke/SubTVec.h>
#include <Kripke/SweepSolver.h>
#include <vector>
#include <stdio.h>


/**
  Run solver iterations.
*/
int Kripke::SteadyStateSolver (Kripke::DataStore &data_store, Grid_Data *grid_data, bool block_jacobi)
{
  KRIPKE_TIMER(data_store, Solve);


  Kripke::Comm const &comm = data_store.getVariable<Kripke::Comm>("comm");
  if(comm.rank() == 0){
    printf("\n");
    printf("Steady State Solve\n");
    printf("==================\n\n");
  }


  // Loop over iterations
  double part_last = 0.0;
  for(int iter = 0;iter < grid_data->niter;++ iter){


    /*
     * Compute the RHS:  rhs = LPlus*S*L*psi + Q
     */

    // Discrete to Moments transformation (phi = L*psi)
    Kripke::Kernel::LTimes(data_store);

    // Compute Scattering Source Term (psi_out = S*phi)
    Kripke::Kernel::scattering(data_store);

    // Compute External Source Term (psi_out = psi_out + Q)
    Kripke::Kernel::source(data_store);

    // Moments to Discrete transformation (rhs = LPlus*psi_out)
    Kripke::Kernel::LPlusTimes(data_store);




    /*
     * Sweep (psi = Hinv*rhs)
     */
    {
      // Create a list of all groups
      std::vector<int> sdom_list(grid_data->subdomains.size());
      for(size_t i = 0;i < grid_data->subdomains.size();++ i){
        sdom_list[i] = i;
      }

      // Sweep everything
      Kripke::SweepSolver(data_store, sdom_list, grid_data, block_jacobi);
    }



    /*
     * Population edit and convergence test
     */
    double part = Kripke::Kernel::population(data_store);
    if(comm.rank() == 0){
      printf("  iter %d: particle count=%e, change=%e\n", iter, part, (part-part_last)/part);
    }
    part_last = part;



  }

  printf("  Solver terminated\n");

  return(0);
}




