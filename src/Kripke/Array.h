/*--------------------------------------------------------------------------
 * Header file for the Grid_Data data structures
 *--------------------------------------------------------------------------*/

#ifndef KRIPKE_ARRAY_H__
#define KRIPKE_ARRAY_H__

struct Array3 {
  Array3();
  ~Array3();

  void resize(int n0, int n1, int n2);
  void dealloc(void);

  inline double **operator[](int idx){
    return data[idx];
  }


  int dim[3];
  std::vector<double> data_linear;
  double ***data;
};

struct Array5 {
  Array5();
  ~Array5();

  void resize(int n0, int n1, int n2, int n3, int n4);

  inline double ****operator[](int idx){
    return data[idx];
  }

  int dim[5];
  std::vector<double> data_linear;
  double *****data;
};

#endif
