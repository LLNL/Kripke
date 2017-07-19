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

#ifndef KRIPKE_SUBDOMAIN_H__
#define KRIPKE_SUBDOMAIN_H__

#include <vector>
#include <Kripke/Layout.h>

// Foreward Decl
struct SubTVec;
struct InputVariables;

namespace Kripke {
  struct QuadraturePoint;
}


/**
 * Contains parameters and variables that describe a single Group Set and
 * Direction Set.
 */
struct Subdomain {
  Subdomain();
  ~Subdomain();

  void setup(int sdom_id, InputVariables *input_vars, int gs, int ds, int zs,
    std::vector<Kripke::QuadraturePoint> &direction_list, Layout *layout);

  // Neighbors
  Neighbor upwind[3];   // Upwind dependencies in x,y,z
  Neighbor downwind[3]; // Downwind neighbors in x,y,z

};

#endif
