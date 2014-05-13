#ifndef KRIPKE_H__
#define KRIPKE_H__

#include<string>
#include<vector>
#include<stdio.h>
#include<cmath>

/* Prototypes */
struct User_Data;

#define KRIPKE_USE_PAPI 1
#define KRIPKE_USE_PERFTOOLS 1

/* driver.c */
void Driver(User_Data *user_data);

/* sweep_solver.c */
int SweepSolver(User_Data *user_data);
int SweepSolver_GroupSet (int group_set, User_Data *user_data);
void CreateBufferInfo(User_Data *user_data);


enum Nesting_Order {
  // Nestings for Psi and Phi
  // D referes to directions OR moments, depending on context
  NEST_GDZ,
  NEST_DGZ,
  NEST_ZDG,
  NEST_DZG,
  NEST_ZGD,
  NEST_GZD,

  // Nestings for L and L+ matrices
  NEST_DNM,
  NEST_NMD
};

inline std::string nestingString(Nesting_Order nesting){
  switch(nesting){
    case NEST_GDZ: return("GDZ");
    case NEST_GZD: return("GZD");
    case NEST_ZDG: return("ZDG");
    case NEST_ZGD: return("ZGD");
    case NEST_DGZ: return("DGZ");
    case NEST_DZG: return("DZG");
    case NEST_DNM: return("DNM");
    case NEST_NMD: return("NMD");
  }
  return("UNKNOWN");
}

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

