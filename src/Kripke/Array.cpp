#include<Kripke/Array.h>
#include<stdio.h>

Array3::Array3() :
  dim({0,0,0}),
  data_linear(0),
  data(NULL)
{
}

Array3::~Array3(){
  dealloc();
}

void Array3::resize(int n0, int n1, int n2){
  dealloc();
  dim[0] = n0;
  dim[1] = n1;
  dim[2] = n2;

  // allocate block of memory
  data_linear.resize(n0*n1*n2);

  // create pointers into that data
  data = new double**[n0];
  double *ptr = &data[0];
  for(int d0 = 0;d0 < dim[0];d0 ++){
    data[d0] = new double*[n1];
    for(int d1 = 0;d1 < dim[1];d1 ++){
      data[d0][d1] = ptr;
      ptr += dim[2];
    }
  }
}


void Array3::dealloc(void){
  if(data == NULL){
    return;
  }

  // deallocate data pointers
  for(int d0 = 0;d0 < dim[0];d0 ++){
    delete[] data[d0];
  }
  delete[] data;

  // clear everything else
  data_linear.clear();
  dim[0] = 0;
  dim[1] = 0;
  dim[2] = 0;
}



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
