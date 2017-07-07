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

#include <Kripke/Initialize.h>

#include <Kripke/Comm.h>
#include <Kripke/Field.h>
#include <Kripke/PartitionSpace.h>
#include <Kripke/Set.h>
#include <Kripke/Timing.h>
#include <Kripke/VarTypes.h>

using namespace Kripke;


static void initializeDecomp(Kripke::DataStore &data_store,
    InputVariables const &input_vars)
{
  // Create a "Comm World"
  Comm *comm = new Comm();
  data_store.addVariable("comm", comm);

  // Create our partitioning
  PartitionSpace *pspace = new PartitionSpace(*comm,
      1,
      1,
      input_vars.npx,
      input_vars.npy,
      input_vars.npz);

  pspace->setup_createSubdomains(
      input_vars.num_groupsets,
      input_vars.num_dirsets,
      input_vars.num_zonesets_dim[0],
      input_vars.num_zonesets_dim[1],
      input_vars.num_zonesets_dim[2]);

  pspace->print();

  data_store.addVariable("pspace", pspace);
}





static void initializeEnergy(Kripke::DataStore &data_store,
    InputVariables const &input_vars)
{

  PartitionSpace &pspace = data_store.getVariable<PartitionSpace>("pspace");

  // Create sets for energy discretization
  size_t ngrp_per_sdom = input_vars.num_groups /
                         pspace.getGlobalNumSubdomains(SPACE_P);

  std::vector<size_t> local_grps(pspace.getNumSubdomains(SPACE_P),
                                 ngrp_per_sdom);

  RangeSet *grp_set = new RangeSet(pspace, SPACE_P, local_grps);

  data_store.addVariable("Set/Group", grp_set);


}





static void initializeDirections(Kripke::DataStore &data_store,
    InputVariables const &input_vars)
{

  PartitionSpace &pspace = data_store.getVariable<PartitionSpace>("pspace");

  // Create sets for angular discretization
  size_t ndir_per_sdom = input_vars.num_directions /
                        pspace.getGlobalNumSubdomains(SPACE_Q);

  std::vector<size_t> local_dirs(pspace.getNumSubdomains(SPACE_Q),
                                 ndir_per_sdom);

  RangeSet *dir_set = new RangeSet(pspace, SPACE_Q, local_dirs);

  data_store.addVariable("Set/Direction", dir_set);


  size_t legendre_order = input_vars.legendre_order;
  size_t num_moments = (legendre_order+1)*(legendre_order+1);
  GlobalRangeSet *moment_set = new GlobalRangeSet(pspace, num_moments);

  data_store.addVariable("Set/Moment", moment_set);

}






static void initializeSpace(Kripke::DataStore &data_store,
    InputVariables const &input_vars)
{
  PartitionSpace &pspace = data_store.getVariable<PartitionSpace>("pspace");



  // Create set for X mesh
  size_t nx_per_sdom = input_vars.nx /
                       pspace.getGlobalNumSubdomains(SPACE_RX);

  std::vector<size_t> local_nx(pspace.getNumSubdomains(SPACE_RX),
                               nx_per_sdom);

  RangeSet *zonex_set = new RangeSet(pspace, SPACE_RX, local_nx);

  data_store.addVariable("Set/ZoneI", zonex_set);




  // Create set for Y mesh
  size_t ny_per_sdom = input_vars.ny /
                       pspace.getGlobalNumSubdomains(SPACE_RY);

  std::vector<size_t> local_ny(pspace.getNumSubdomains(SPACE_RY),
                               ny_per_sdom);

  RangeSet *zoney_set = new RangeSet(pspace, SPACE_RY, local_ny);

  data_store.addVariable("Set/ZoneJ", zoney_set);




  // Create set for Z mesh
  size_t nz_per_sdom = input_vars.nz /
                       pspace.getGlobalNumSubdomains(SPACE_RZ);

  std::vector<size_t> local_nz(pspace.getNumSubdomains(SPACE_RZ),
                               nz_per_sdom);

  RangeSet *zonez_set = new RangeSet(pspace, SPACE_RZ, local_nz);

  data_store.addVariable("Set/ZoneK", zonez_set);



  // Create a total set of zones in the problem
  ProductSet<3> *zone_set = new ProductSet<3>(pspace, SPACE_R,
      *zonex_set, *zoney_set, *zonez_set);

  data_store.addVariable("Set/Zone", zone_set);
}






static void initializeData(Kripke::DataStore &data_store,
    InputVariables const &)
{

  PartitionSpace &pspace = data_store.getVariable<PartitionSpace>("pspace");


  // Create a set to span angular the flux
  Set const &dir_set   = data_store.getVariable<Set>("Set/Direction");
  Set const &group_set = data_store.getVariable<Set>("Set/Group");
  Set const &zone_set  = data_store.getVariable<Set>("Set/Zone");
  ProductSet<3> *flux_set = new ProductSet<3>(pspace, SPACE_PQR,
      dir_set, group_set, zone_set);

  data_store.addVariable("Set/Flux", flux_set);

  // Create Solution and RHS fields
  data_store.addVariable("psi", new Field_Flux(*flux_set));
  data_store.addVariable("rhs", new Field_Flux(*flux_set));




  // Create a set to span moments of the angular flux
  Set const &moment_set   = data_store.getVariable<Set>("Set/Moment");
  ProductSet<3> *fluxmoment_set = new ProductSet<3>(pspace, SPACE_PR,
        moment_set, group_set, zone_set);

  data_store.addVariable("Set/FluxMoment", fluxmoment_set);


  // Create flux moment and source moment fields
  data_store.addVariable("phi",     new Field_Flux(*fluxmoment_set));
  data_store.addVariable("phi_out", new Field_Flux(*fluxmoment_set));



  Comm default_comm;
  if(default_comm.rank() == 0){

    unsigned long flux_size = flux_set->globalSize();
    unsigned long moment_size = fluxmoment_set->globalSize();
    unsigned long total_size = (2*moment_size + 2*flux_size);

    printf("\n");
    printf("  Variable     Num Elements      Megabytes\n");
    printf("  --------     ------------      ---------\n");
    printf("  psi          %12lu   %12.3lf\n",
        flux_size,
        (double)flux_size*8.0/1024.0/1024.0);
    printf("  rhs          %12lu   %12.3lf\n",
        flux_size,
        (double)flux_size*8.0/1024.0/1024.0);

    printf("  phi          %12lu   %12.3lf\n",
        moment_size,
        (double)moment_size*8.0/1024.0/1024.0);
    printf("  phi_out      %12lu   %12.3lf\n",
        moment_size,
        (double)moment_size*8.0/1024.0/1024.0);

    printf("  TOTAL        %12lu   %12.3lf\n",
        total_size,
        (double)total_size*8.0/1024.0/1024.0);
  }

}





void Kripke::initializeDataStore(Kripke::DataStore &data_store,
    InputVariables const &input_vars)
{

  Comm default_comm;

  if(default_comm.rank() == 0){
    printf("\nInitializing\n");
    printf("============\n\n");
  }

  // Create and start a timing object
  data_store.addVariable("timing", new Kripke::Timing());
  KRIPKE_TIMER(data_store, Initialize);


  // Create parallel and subdomain decomposition
  initializeDecomp(data_store, input_vars);

  // Create energy discretization
  initializeEnergy(data_store, input_vars);

  // Create angular discretization, quadrature set and L/L+ matrices
  initializeDirections(data_store, input_vars);

  // Create a spatial mesh, and paint it with materials
  initializeSpace(data_store, input_vars);

  // Create cross sections and transfer matrix
  initializeData(data_store, input_vars);
}
