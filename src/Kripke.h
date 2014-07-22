/**
 * This is the main header file for the Kripke Mini-App.
 */
#ifndef KRIPKE_H__
#define KRIPKE_H__

#include<string>
#include<vector>
#include<stdio.h>
#include<cmath>
#include<strings.h>

// Make sure that there's openmp support, otherwise error out
#if KRIPKE_USE_OPENMP
#ifndef _OPENMP
#error "OpenMP selected for build, but OpenMP is not available"
#endif
#endif


// Forward Decl
struct User_Data;

#define KRESTRICT __restrict__


// In Kripke/Sweep_Solver.cpp
int SweepSolver(User_Data *user_data);
int SweepSolver_GroupSet (int group_set, User_Data *user_data);
void CreateBufferInfo(User_Data *user_data);

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
  NEST_ZGD,

  // Nestings for L and L+ matrices
  NEST_DNM,
  NEST_NMD
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

    case NEST_DNM: return("DNM");
    case NEST_NMD: return("NMD");
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


/**
 * Compares two vectors for differences.
 * Used in testing suite.
 */
inline bool compareVector(std::string const &name,
    std::vector<double> const &a,
    std::vector<double> const &b, double tol, bool verbose){

  if(a.size() != b.size()){
    if(verbose){
      printf("Vectors are different lengths: %ld, %ld\n",
          (long)a.size(), (long)b.size());
    }
    return true;
  }

  bool is_diff = false;
  for(size_t i = 0;i < a.size();++i){
    if(std::abs(a[i]-b[i]) > tol){
      is_diff = true;
      if(verbose){
        printf("%s[%d]:%e, %e [%e]\n",
            name.c_str(), (int)i,
            a[i], b[i], std::abs(a[i]-b[i]));
        is_diff = true;
      }
      else{
        break;
      }
    }
  }

  return is_diff;
}

/**
 * Compares two scalars for differences.
 * Used in testing suite.
 */
inline bool compareScalar(std::string const &name,
    double a, double b, double tol, bool verbose){

  if(std::abs(a-b) > tol){
    if(verbose){
      printf("%s:%e, %e [%e]\n",
          name.c_str(),
          a, b, std::abs(a-b));
    }
    return true;
  }
  return false;
}

#endif

