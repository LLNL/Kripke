#ifndef KRIPKE_SUBTVEC_H__
#define KRIPKE_SUBTVEC_H__

#include <Kripke/Kernel.h>
#include <algorithm>
#include <vector>
#include <stdlib.h>

/* A transport vector (used for Psi and Phi, RHS, etc.)
 *
 *It is templated to re-order the nestings
 *
 * 0: outer nesting
 * 1: middle nesting
 * 2: inner nesting
 *
 */
struct SubTVec {
  SubTVec(Nesting_Order nesting, int ngrps, int ndir_mom, int nzones):
    groups(ngrps),
    directions(ndir_mom),
    zones(nzones),
    elements(groups*directions*zones),
    data(NULL),
    data_linear(elements)
  {
    // setup nesting order
    int int_to_ext[3];
    switch(nesting){
      case NEST_GDZ:
        int_to_ext[0] = 0;
        int_to_ext[1] = 1;
        int_to_ext[2] = 2;
        break;
      case NEST_GZD:
        int_to_ext[0] = 0;
        int_to_ext[2] = 1;
        int_to_ext[1] = 2;
        break;
      case NEST_DZG:
        int_to_ext[1] = 0;
        int_to_ext[2] = 1;
        int_to_ext[0] = 2;
        break;
      case NEST_DGZ:
        int_to_ext[1] = 0;
        int_to_ext[0] = 1;
        int_to_ext[2] = 2;
        break;
      case NEST_ZDG:
        int_to_ext[2] = 0;
        int_to_ext[1] = 1;
        int_to_ext[0] = 2;
        break;
      case NEST_ZGD:
        int_to_ext[2] = 0;
        int_to_ext[0] = 1;
        int_to_ext[1] = 2;
        break;
    }

    // setup dimensionality
    int size_ext[3];
    size_ext[0] = groups;
    size_ext[1] = directions;
    size_ext[2] = zones;

    // map to internal indices
    for(int i = 0; i < 3; ++i){
      ext_to_int[i] = int_to_ext[i];
    }
    for(int i = 0; i < 3; ++i){
      size_int[ext_to_int[i]] = size_ext[i];
    }

    data = new double**[size_int[0]];
    for(int a = 0; a < size_int[0]; ++a){
      data[a] = new double*[size_int[1]];
      for(int b = 0; b < size_int[1]; ++b){
        data[a][b] = &data_linear[0] + a * size_int[1]*size_int[2] + b *
                     size_int[2];
      }
    }
  }


  /** ALIASING version */
  SubTVec(Nesting_Order nesting, int ngrps, int ndir_mom, int nzones, double *ptr):
    groups(ngrps),
    directions(ndir_mom),
    zones(nzones),
    elements(groups*directions*zones),
    data(NULL),
    data_linear(0)
  {
    // setup nesting order
    int int_to_ext[3];
    switch(nesting){
      case NEST_GDZ:
        int_to_ext[0] = 0;
        int_to_ext[1] = 1;
        int_to_ext[2] = 2;
        break;
      case NEST_GZD:
        int_to_ext[0] = 0;
        int_to_ext[2] = 1;
        int_to_ext[1] = 2;
        break;
      case NEST_DZG:
        int_to_ext[1] = 0;
        int_to_ext[2] = 1;
        int_to_ext[0] = 2;
        break;
      case NEST_DGZ:
        int_to_ext[1] = 0;
        int_to_ext[0] = 1;
        int_to_ext[2] = 2;
        break;
      case NEST_ZDG:
        int_to_ext[2] = 0;
        int_to_ext[1] = 1;
        int_to_ext[0] = 2;
        break;
      case NEST_ZGD:
        int_to_ext[2] = 0;
        int_to_ext[0] = 1;
        int_to_ext[1] = 2;
        break;
    }

    // setup dimensionality
    int size_ext[3];
    size_ext[0] = groups;
    size_ext[1] = directions;
    size_ext[2] = zones;

    // map to internal indices
    for(int i = 0; i < 3; ++i){
      ext_to_int[i] = int_to_ext[i];
    }
    for(int i = 0; i < 3; ++i){
      size_int[ext_to_int[i]] = size_ext[i];
    }

    data = new double**[size_int[0]];
    for(int a = 0; a < size_int[0]; ++a){
      data[a] = new double*[size_int[1]];
      for(int b = 0; b < size_int[1]; ++b){
        data[a][b] = ptr + a * size_int[1]*size_int[2] + b *
                     size_int[2];
      }
    }
  }

  ~SubTVec(){
    for(int a = 0; a < size_int[0]; ++a){
      delete[] data[a];
    }
    delete[] data;
  }

  // These are NOT efficient.. just used to re-stride data for comparisons
  inline double &operator()(int g, int d, int z) {
    int idx[3];
    idx[ext_to_int[0]] = g;
    idx[ext_to_int[1]] = d;
    idx[ext_to_int[2]] = z;
    return(data[idx[0]][idx[1]][idx[2]]);
  }
  inline double operator()(int g, int d, int z) const {
    int idx[3];
    idx[ext_to_int[0]] = g;
    idx[ext_to_int[1]] = d;
    idx[ext_to_int[2]] = z;
    return(data[idx[0]][idx[1]][idx[2]]);
  }

  inline double sum(void) const {
    double s = 0.0;
    for(size_t i = 0;i < data_linear.size();++ i){
      s+= data_linear[i];
    }
    return s;
  }

  inline void clear(double v){
    std::fill(data_linear.begin(), data_linear.end(), v);
  }

  inline void randomizeData(void){
    for(int i = 0;i < data_linear.size();++ i){
      data_linear[i] = drand48();
    }
  }

  inline void copy(SubTVec const &b){
    for(int g = 0;g < groups;++ g){
      for(int d = 0;d < directions; ++ d){
        for(int z = 0;z < zones;++ z){
          // Copy using abstract indexing
          (*this)(g,d,z) = b(g,d,z);
        }
      }
    }
  }

  inline bool compare(std::string const &name, SubTVec const &b,
      double tol, bool verbose){

    bool is_diff = false;
    for(int g = 0;g < groups;++ g){
      for(int d = 0;d < directions; ++ d){
        for(int z = 0;z < zones;++ z){
          // Copy using abstract indexing
          double err = std::abs((*this)(g,d,z) - b(g,d,z));
          if(err > tol){
            is_diff = true;
            if(verbose){
              printf("%s[g=%d, d=%d, z=%d]: |%e - %e| = %e\n",
                  name.c_str(), g,d,z, (*this)(g,d,z), b(g,d,z), err);
            }
          }
        }
      }
    }
    return is_diff;
  }

  int ext_to_int[3]; // external index to internal index mapping
  int size_int[3]; // size of each dimension in internal indices

  int groups, directions, zones, elements;
  double ***data;
  std::vector<double> data_linear;
};


#endif
