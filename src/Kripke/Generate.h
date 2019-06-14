//
// Copyright (c) 2014-19, Lawrence Livermore National Security, LLC
// and Kripke project contributors. See the COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
//

#ifndef KRIPKE_GENERATE_H__
#define KRIPKE_GENERATE_H__

#include <Kripke.h>
#include <Kripke/Core/DataStore.h>
#include <Kripke/InputVariables.h>

namespace Kripke {

  /**
   * Takes an Input_Variables object and generates a problem in the DataStore
   */
  void generateProblem(Kripke::Core::DataStore &data_store,
      InputVariables const &input_variables);
  

  namespace Generate {


    void generateDecomp(Kripke::Core::DataStore &data_store,
        InputVariables const &input_variables);

    void generateEnergy(Kripke::Core::DataStore &data_store,
          InputVariables const &input_variables);

    void generateQuadrature(Kripke::Core::DataStore &data_store,
          InputVariables const &input_variables);

    void generateSpace(Kripke::Core::DataStore &data_store,
          InputVariables const &input_variables);

    void generateData(Kripke::Core::DataStore &data_store,
          InputVariables const &input_variables);

  }

}

#endif
