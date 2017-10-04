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

#include <Kripke/Generate.h>

#include <Kripke/Core/Comm.h>
#include <Kripke/Core/Field.h>
#include <Kripke/Kernel.h>
#include <Kripke/Core/PartitionSpace.h>
#include <Kripke/Core/Set.h>
#include <Kripke/Timing.h>
#include <Kripke/VarTypes.h>

using namespace Kripke;
using namespace Kripke::Core;


void Kripke::Generate::generateData(Kripke::Core::DataStore &data_store,
    InputVariables const &input_vars)
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
  data_store.addVariable("phi",     new Field_Moments(*fluxmoment_set));
  data_store.addVariable("phi_out", new Field_Moments(*fluxmoment_set));


  // Create "plane data" to hold face-centered values while sweeping
  Set const &zonei_set = data_store.getVariable<Set>("Set/ZoneI");
  Set const &zonej_set = data_store.getVariable<Set>("Set/ZoneJ");
  Set const &zonek_set = data_store.getVariable<Set>("Set/ZoneK");
  Set const &iplane_set = data_store.newVariable<ProductSet<4>>("Set/IPlane", pspace, SPACE_PQR, dir_set, group_set, zonej_set, zonek_set);
  Set const &jplane_set = data_store.newVariable<ProductSet<4>>("Set/JPlane", pspace, SPACE_PQR, dir_set, group_set, zonei_set, zonek_set);
  Set const &kplane_set = data_store.newVariable<ProductSet<4>>("Set/KPlane", pspace, SPACE_PQR, dir_set, group_set, zonei_set, zonej_set);
  data_store.newVariable<Field_IPlane>("i_plane", iplane_set);
  data_store.newVariable<Field_JPlane>("j_plane", jplane_set);
  data_store.newVariable<Field_KPlane>("k_plane", kplane_set);

  // Create a set to span scattering transfer matrix
  Set const &material_set   = data_store.getVariable<Set>("Set/Material");
  Set const &legendre_set   = data_store.getVariable<Set>("Set/Legendre");
  Set const &global_group_set = data_store.getVariable<Set>("Set/GlobalGroup");
  ProductSet<4> *sigs_set = new ProductSet<4>(pspace, SPACE_NULL,
      material_set, legendre_set, global_group_set, global_group_set);

  data_store.addVariable("Set/SigmaS", sigs_set);


  // Create storage for the scattering transfer matrix
  data_store.addVariable("data/sigs", new Field_SigmaS(*sigs_set));
  auto &field_sigs = data_store.getVariable<Field_SigmaS>("data/sigs");

  // Assign basic diagonal data to matrix
  for(auto sdom_id : field_sigs.getWorkList()){

    // Zero out entire matrix
    auto sigs_ptr = field_sigs.getData(sdom_id);
    size_t sigs_size = field_sigs.size(sdom_id);
    for(size_t i = 0;i < sigs_size;++ i){
      sigs_ptr[i] = 0.0;
    }

    // Assign diagonal to the user input for each material
    // Assume each group has same behavior
    auto sigs = field_sigs.getView(sdom_id);
    int global_num_groups = global_group_set.size(sdom_id);
    Legendre n{0};
    for(Material mat{0};mat < 3;++ mat){
      for(GlobalGroup g{0};g < global_num_groups;++ g){
        sigs(mat, n, g, g) = input_vars.sigs[*mat];
      }
    }
  }



}




