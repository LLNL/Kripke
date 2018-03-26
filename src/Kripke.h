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

#ifndef KRIPKE_H__
#define KRIPKE_H__

#include <KripkeConfig.h>

#include <RAJA/RAJA.hpp>

#include <string>
#include <vector>
#include <stdio.h>
#include <cmath>
#include <strings.h>
#include <stdlib.h>

// Make sure that there's openmp support, otherwise error out
#ifdef KRIPKE_USE_OPENMP
#ifndef _OPENMP
#error "OpenMP selected for build, but OpenMP is not available"
#endif
#endif

#ifdef KRIPKE_USE_MPI
#include <mpi.h>
#endif

// Forward Decl
struct Grid_Data;

#define KRESTRICT __restrict__

#ifdef KRIPKE_USE_MPI
#define KRIPKE_ABORT(...) \
  printf(__VA_ARGS__); \
  MPI_Abort(MPI_COMM_WORLD, 1);
#else
#define KRIPKE_ABORT(...) \
  printf(__VA_ARGS__); \
  exit(1);
#endif


#define KRIPKE_ASSERT(EXPR, ...) \
  if(!(EXPR)){\
    KRIPKE_ABORT("Assertion Failed: " __VA_ARGS__); \
  }
   

#define KRIPKE_LAMBDA [=] RAJA_HOST_DEVICE

namespace Kripke {

  /**
   * Index used to specify a local subdomain
   */
  RAJA_INDEX_VALUE(SdomId, "SdomId");


  /**
   * Index used to specify a global subdomain
   */
  RAJA_INDEX_VALUE(GlobalSdomId, "GlobalSdomId");


}



/**
 * Tags for choosing which data nesting to be chosen
 */
enum Nesting_Order {
  // Nestings for Psi and Phi
  // D referes to directions OR moments, depending on context
  NEST_DGZ,
  NEST_DZG,
  NEST_GDZ,
  NEST_GZD,
  NEST_ZDG,
  NEST_ZGD
};


/**
  Tags for which parallel algorithm to use.
*/
enum ParallelMethod {
  PMETHOD_SWEEP,
  PMETHOD_BJ
};

/**
 * Converts a nesting tag to a human-readable string.
 */
inline std::string nestingString(Nesting_Order nesting){
  switch(nesting){
    case NEST_DGZ: return("DGZ");
    case NEST_DZG: return("DZG");
    case NEST_GDZ: return("GDZ");
    case NEST_GZD: return("GZD");
    case NEST_ZDG: return("ZDG");
    case NEST_ZGD: return("ZGD");
  }
  return("UNKNOWN");
}

/**
 * Converts a string (eg. from command line) to a nesting tag.
 */
inline Nesting_Order nestingFromString(std::string const &str){
  for(int i = 0;i < 6;++ i){
    if(!strcasecmp(str.c_str(), nestingString((Nesting_Order)i).c_str())){
      return (Nesting_Order)i;
  }
 }
  return (Nesting_Order)-1;
}




#endif

