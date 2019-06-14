//
// Copyright (c) 2014-19, Lawrence Livermore National Security, LLC
// and Kripke project contributors. See the COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
//

#include <Kripke/Generate.h>

#include <Kripke/Core/Comm.h>
#include <Kripke/Core/Field.h>
#include <Kripke/Kernel.h>
#include <Kripke/Core/PartitionSpace.h>
#include <Kripke/Core/Set.h>
#include <Kripke/Timing.h>
#include <Kripke/VarTypes.h>
#include <string>
#include <vector>

using namespace Kripke;
using namespace Kripke::Core;


void Kripke::generateProblem(Kripke::Core::DataStore &data_store,
    InputVariables const &input_vars)
{

  Comm default_comm;

  if(default_comm.rank() == 0){
    printf("\nGenerating Problem\n");
    printf("==================\n\n");
  }

  // Create and start a timing object
  data_store.addVariable("timing", new Kripke::Timing());
  KRIPKE_TIMER(data_store, Generate);


  // Create parallel and subdomain decomposition
  Generate::generateDecomp(data_store, input_vars);

  // Create energy discretization
  Generate::generateEnergy(data_store, input_vars);

  // Create angular discretization, quadrature set and L/L+ matrices
  Generate::generateQuadrature(data_store, input_vars);

  // Create a spatial mesh, and paint it with materials
  Generate::generateSpace(data_store, input_vars);

  // Create cross sections and transfer matrix
  Generate::generateData(data_store, input_vars);



  // Display all of the fields that were created, and what their sizes are
  if(default_comm.rank() == 0){

    // Collect variables that are Fields of doubles
    std::vector<std::string> field_names;
    for(auto const &var_name : data_store.getVariableList()){
      if(data_store.isVariableType<FieldStorage<double>>(var_name)){
        field_names.push_back(var_name);
      }
    }
    std::sort(field_names.begin(), field_names.end());

    printf("\n");
    printf("  Memory breakdown of Field variables:\n");
    printf("  Field Variable            Num Elements    Megabytes\n");
    printf("  --------------            ------------    ---------\n");

    unsigned long total_size = 0;
    for(auto const &field_name : field_names){

      unsigned long field_size = data_store.getVariable<FieldStorage<double>>(field_name).getSet().globalSize();
      total_size += field_size;

      printf("  %-24s  %12lu %12.3lf\n",
          field_name.c_str(),
          field_size,
          (double)field_size*8.0/1024.0/1024.0);
    }

    printf("  --------                  ------------    ---------\n");
    printf("  TOTAL                     %12lu %12.3lf\n",
        total_size,
        (double)total_size*8.0/1024.0/1024.0);

    printf("\n");
    printf("  Generation Complete!\n");
  }
}
