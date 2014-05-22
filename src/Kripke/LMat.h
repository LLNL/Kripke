#ifndef KRIPKE_LMAT_H__
#define KRIPKE_LMAT_H__

#include <Kripke/Kernel.h>
#include <algorithm>
#include <vector>
#include <stdlib.h>

/*
 * An L or L+ matrix (used for computing moments)
 */
struct LMat {
  LMat(Nesting_Order nesting, int dims, int moments, int dirs):
    num_dims(dims),
    num_directions(dirs),
    num_moments(moments),
    data(NULL),
    data_linear()
  {
    // setup nesting order
    int int_to_ext[3];
    switch(nesting){
      case NEST_NMD:
        int_to_ext[0] = 0;
        int_to_ext[1] = 1;
        int_to_ext[2] = 2;
        break;
      case NEST_DNM:
        int_to_ext[2] = 0;
        int_to_ext[0] = 1;
        int_to_ext[1] = 2;
        break;
      default:
        throw "Invalid nesting for L or L+ matrix";
    }

    // Figure out how many total M's there are
    int largest_m;
    switch(dims){
      case 1:
        largest_m = 1;
        total_moments = moments;
        break;
      case 2:
        largest_m = moments+1;
        total_moments = (largest_m+1)*largest_m / 2;
        break;
      case 3:
        largest_m = moments*2+1;
        total_moments = moments*moments;
        break;
      default:
        throw "Invalid number of dimensions";
    }

    // Allocate storage
    size_t elements = total_moments * dirs;
    data_linear.resize(elements, 0.0);

    // setup dimensionality
    int size_ext[3];
    size_ext[0] = num_moments;
    size_ext[1] = largest_m;
    size_ext[2] = num_directions;

    // map to internal indices
    for(int i = 0; i < 3; ++i){
      ext_to_int[i] = int_to_ext[i];
    }
    for(int i = 0; i < 3; ++i){
      size_int[ext_to_int[i]] = size_ext[i];
    }

    // Setup Nested Data

    // Nested NMD
    if(nesting == NEST_NMD){
      data = new double**[num_moments];
      double *ptr = &data_linear[0];
      for(int n = 0; n < num_moments; ++n){
        int num_m = numM(n);
        data[n] = new double*[num_m];

        for(int m = 0;m < num_m;++ m){
          data[n][m] = ptr;
          ptr += num_directions;
        }
      }
    }
    else{ // Nested DNM
      data = new double**[num_directions];
      double *ptr = &data_linear[0];
      for(int d = 0;d < num_directions;++ d){
        data[d] = new double*[num_moments];
        for(int n = 0; n < num_moments; ++n){
          data[d][n] = ptr;
          ptr += numM(n);
        }
      }
    }
  }

  ~LMat(){
    for(int a = 0; a < size_int[0]; ++a){
      delete[] data[a];
    }
    delete[] data;
  }

  inline int numM(int n) const {
    switch(num_dims){
      case 1:
        return 1;
      case 2:
        return 1+n;
      case 3:
        return 1+2*n;
    }
    return 0;
  }

  inline double* ptr(void){
    return &data_linear[0];
  }

  // These are NOT efficient.. just used to re-stride data for comparisons
  inline double &operator()(int n, int m, int d) {
    int idx[3];
    idx[ext_to_int[0]] = n;
    idx[ext_to_int[1]] = m;
    idx[ext_to_int[2]] = d;
    return(data[idx[0]][idx[1]][idx[2]]);
  }
  inline double operator()(int n, int m, int d) const {
    int idx[3];
    idx[ext_to_int[0]] = n;
    idx[ext_to_int[1]] = m;
    idx[ext_to_int[2]] = d;
    return(data[idx[0]][idx[1]][idx[2]]);
  }

  inline void clear(double v){
    std::fill(data_linear.begin(), data_linear.end(), v);
  }

  inline void randomizeData(void){
    for(int i = 0;i < data_linear.size();++ i){
      data_linear[i] = drand48();
    }
  }

  inline void copy(LMat const &b){
    for(int n = 0;n < num_moments;++ n){
      int num_m = numM(n);
      for(int m = 0;m < num_m;++ m){
        for(int d = 0;d < num_directions;++ d){
          (*this)(n,m,d) = b(n,m,d);
        }
      }
    }
  }

  inline bool compare(std::string const &name, LMat const &b,
      double tol, bool verbose){
    bool is_diff = false;
    for(int n = 0;n < num_moments;++ n){
      int num_m = numM(n);
      for(int m = 0;m < num_m;++ m){
        for(int d = 0;d < num_directions;++ d){
          // Copy using abstract indexing
          double err = std::abs((*this)(n,m,d) - b(n,m,d));
          if(err > tol){
            is_diff = true;
            if(verbose){
              printf("%s[n=%d, m=%d, d=%d]: |%e - %e| = %e\n",
                  name.c_str(), n,m,d, (*this)(n,m,d), b(n,m,d), err);
            }
          }
        }
      }
    }
    return is_diff;
  }

  int ext_to_int[3]; // external index to internal index mapping
  int size_int[3]; // size of each dimension in internal indices

  int num_dims, num_moments, num_directions, total_moments;
  double ***data;
  std::vector<double> data_linear;
};


#endif
