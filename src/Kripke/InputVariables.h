//
// Copyright (c) 2014-19, Lawrence Livermore National Security, LLC
// and Kripke project contributors. See the COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
//

#ifndef KRIPKE_INPUT_VARIABLES_H__
#define KRIPKE_INPUT_VARIABLES_H__

#include <Kripke.h>
#include <Kripke/ArchLayout.h>

/**
 * This structure defines the input parameters to setup a problem.
 */

struct InputVariables {
  InputVariables();
  
  bool checkValues(void) const;
  
  // Problem Description
  int nx, ny, nz;               // Number of spatial zones in x,y,z
  int num_directions;           // Total number of directions
  int num_groups;               // Total number of energy groups
  int legendre_order;           // Scattering order (number Legendre coeff's - 1)
  int quad_num_polar;           // Number of polar quadrature points (0 for dummy)
  int quad_num_azimuthal;       // Number of azimuthal quadrature points (0 for dummy)

  // On-Node Options
  Kripke::ArchLayoutV al_v;     // Data layout and architecture selection
  
  // Parallel Decomp
  int npx, npy, npz;            // The number of processors in x,y,z
  int num_dirsets;              // Number of direction sets
  int num_groupsets;            // Number of energy group sets
  int num_zonesets_dim[3];      // Number of zoneset in x, y, z  
  
  // Physics and Solver Options
  int niter;                    // number of solver iterations to run
  ParallelMethod parallel_method;
  double sigt[3];               // total cross section for 3 materials
  double sigs[3];               // total scattering cross section for 3 materials
  int num_material_subsamples;  // number of subsamples in each dimension for mesh painting
  
  // Output Options
  std::string run_name;         // Name to use when generating output files
};

#endif
