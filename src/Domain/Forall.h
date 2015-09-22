#ifndef DOMAIN_FORALL_H__
#define DOMAIN_FORALL_H__

struct seq_pol{};
struct omp_pol{};

template<typename BODY>
inline void forall(seq_pol const &pol, int start, int end, BODY const &body){
  for(int i = start;i < end;++ i){
    body(i);
  }
}

template<typename BODY>
inline void forall(omp_pol const &pol, int start, int end, BODY const &body){
#ifdef KRIPKE_USE_OPENMP
#pragma omp parallel for
#endif
  for(int i = start;i < end;++ i){
    body(i);
  }
}

// TODO: How does Kokkos deal with capture by value and constness of Views?
template<typename POL, typename BODY>
inline void forall(int start, int end, BODY const &body){
  forall(POL(), start, end, body);
} 



template<typename POL, typename BODY>
inline void forall2(LAYOUT_IJ const, int end_i, int end_j, BODY const &body){  
  forall<typename POL::pol_i>(0, end_i, [&](int i){
    forall<typename POL::pol_j>(0, end_j, [&](int j){
      body(i,j);
    });
  });
}

template<typename POL, typename BODY>
inline void forall2(LAYOUT_JI const, int end_i, int end_j, BODY const &body){  
  forall<typename POL::pol_j>(0, end_j, [&](int j){
    forall<typename POL::pol_i>(0, end_i, [&](int i){
      body(i,j);
    });
  });
}

template<typename POL, typename BODY>
inline void forall2(int end_i, int end_j, BODY const &body){
  typedef typename POL::layout L;
  forall2<POL, BODY>(L(), end_i, end_j, body);
}



template<typename POL, typename BODY>
inline void forall4(LAYOUT_IJKL const, int end_i, int end_j, int end_k, int end_l, BODY const &body){
  forall<typename POL::pol_i>(0, end_i, [&](int i){
    forall<typename POL::pol_j>(0, end_j, [&](int j){
      forall<typename POL::pol_k>(0, end_k, [&](int k){
        forall<typename POL::pol_l>(0, end_l, [&](int l){
          body(i,j,k,l);  
        });
      });
    });
  });
}

template<typename POL, typename BODY>
inline void forall4(LAYOUT_JIKL const, int end_i, int end_j, int end_k, int end_l, BODY const &body){  
  forall<typename POL::pol_j>(0, end_j, [&](int j){
    forall<typename POL::pol_i>(0, end_i, [&](int i){
      forall<typename POL::pol_k>(0, end_k, [&](int k){
        forall<typename POL::pol_l>(0, end_l, [&](int l){
          body(i,j,k,l);  
        });
      });
    });
  });
}

template<typename POL, typename BODY>
inline void forall4(LAYOUT_JKIL const, int end_i, int end_j, int end_k, int end_l, BODY const &body){  
  forall<typename POL::pol_j>(0, end_j, [&](int j){    
    forall<typename POL::pol_k>(0, end_k, [&](int k){
      forall<typename POL::pol_i>(0, end_i, [&](int i){
        forall<typename POL::pol_l>(0, end_l, [&](int l){
          body(i,j,k,l);  
        });
      });
    });
  });
}

template<typename POL, typename BODY>
inline void forall4(LAYOUT_JKLI const, int end_i, int end_j, int end_k, int end_l, BODY const &body){  
  forall<typename POL::pol_j>(0, end_j, [&](int j){    
    forall<typename POL::pol_k>(0, end_k, [&](int k){
      forall<typename POL::pol_l>(0, end_l, [&](int l){
        forall<typename POL::pol_i>(0, end_i, [&](int i){
          body(i,j,k,l);  
        });
      });
    });
  });
}

template<typename POL, typename BODY>
inline void forall4(LAYOUT_IJLK const, int end_i, int end_j, int end_k, int end_l, BODY const &body){
  forall<typename POL::pol_i>(0, end_i, [&](int i){
    forall<typename POL::pol_j>(0, end_j, [&](int j){          
      forall<typename POL::pol_l>(0, end_l, [&](int l){
        forall<typename POL::pol_k>(0, end_k, [&](int k){
          body(i,j,k,l);  
        });
      });
    });
  });
}

template<typename POL, typename BODY>
inline void forall4(LAYOUT_ILJK const, int end_i, int end_j, int end_k, int end_l, BODY const &body){
  forall<typename POL::pol_i>(0, end_i, [&](int i){
    forall<typename POL::pol_l>(0, end_l, [&](int l){
      forall<typename POL::pol_j>(0, end_j, [&](int j){          
        forall<typename POL::pol_k>(0, end_k, [&](int k){
          body(i,j,k,l);  
        });
      });
    });
  });
}

template<typename POL, typename BODY>
inline void forall4(LAYOUT_JILK const, int end_i, int end_j, int end_k, int end_l, BODY const &body){
  forall<typename POL::pol_j>(0, end_j, [&](int j){          
    forall<typename POL::pol_i>(0, end_i, [&](int i){
      forall<typename POL::pol_l>(0, end_l, [&](int l){
        forall<typename POL::pol_k>(0, end_k, [&](int k){
          body(i,j,k,l);  
        });
      });
    });
  });
}

template<typename POL, typename BODY>
inline void forall4(LAYOUT_KIJL const, int end_i, int end_j, int end_k, int end_l, BODY const &body){
  forall<typename POL::pol_k>(0, end_k, [&](int k){
    forall<typename POL::pol_i>(0, end_i, [&](int i){
      forall<typename POL::pol_j>(0, end_j, [&](int j){          
        forall<typename POL::pol_l>(0, end_l, [&](int l){        
          body(i,j,k,l);  
        });
      });
    });
  });
}

template<typename POL, typename BODY>
inline void forall4(LAYOUT_KJIL const, int end_i, int end_j, int end_k, int end_l, BODY const &body){
  forall<typename POL::pol_k>(0, end_k, [&](int k){    
    forall<typename POL::pol_j>(0, end_j, [&](int j){          
      forall<typename POL::pol_i>(0, end_i, [&](int i){
        forall<typename POL::pol_l>(0, end_l, [&](int l){        
          body(i,j,k,l);  
        });
      });
    });
  });
}


template<typename POL, typename BODY>
inline void forall4(LAYOUT_KLIJ const, int end_i, int end_j, int end_k, int end_l, BODY const &body){
  forall<typename POL::pol_k>(0, end_k, [&](int k){
    forall<typename POL::pol_l>(0, end_l, [&](int l){        
      forall<typename POL::pol_i>(0, end_i, [&](int i){
        forall<typename POL::pol_j>(0, end_j, [&](int j){                  
          body(i,j,k,l);  
        });
      });
    });
  });
}

template<typename POL, typename BODY>
inline void forall4(LAYOUT_KLJI const, int end_i, int end_j, int end_k, int end_l, BODY const &body){
  forall<typename POL::pol_k>(0, end_k, [&](int k){
    forall<typename POL::pol_l>(0, end_l, [&](int l){              
      forall<typename POL::pol_j>(0, end_j, [&](int j){                  
        forall<typename POL::pol_i>(0, end_i, [&](int i){
          body(i,j,k,l);  
        });
      });
    });
  });
}


template<typename POL, typename BODY>
inline void forall4(LAYOUT_LIJK const, int end_i, int end_j, int end_k, int end_l, BODY const &body){  
  forall<typename POL::pol_l>(0, end_l, [&](int l){        
    forall<typename POL::pol_i>(0, end_i, [&](int i){
      forall<typename POL::pol_j>(0, end_j, [&](int j){
        forall<typename POL::pol_k>(0, end_k, [&](int k){
          body(i,j,k,l);  
        });
      });
    });
  });
}

template<typename POL, typename BODY>
inline void forall4(LAYOUT_LJIK const, int end_i, int end_j, int end_k, int end_l, BODY const &body){  
  forall<typename POL::pol_l>(0, end_l, [&](int l){           
    forall<typename POL::pol_j>(0, end_j, [&](int j){
      forall<typename POL::pol_i>(0, end_i, [&](int i){
        forall<typename POL::pol_k>(0, end_k, [&](int k){
          body(i,j,k,l);  
        });
      });
    });
  });
}


template<typename POL, typename BODY>
inline void forall4(LAYOUT_LJKI const, int end_i, int end_j, int end_k, int end_l, BODY const &body){  
  forall<typename POL::pol_l>(0, end_l, [&](int l){           
    forall<typename POL::pol_j>(0, end_j, [&](int j){     
      forall<typename POL::pol_k>(0, end_k, [&](int k){
        forall<typename POL::pol_i>(0, end_i, [&](int i){
          body(i,j,k,l);  
        });
      });
    });
  });
}


template<typename POL, typename BODY>
inline void forall4(LAYOUT_LKIJ const, int end_i, int end_j, int end_k, int end_l, BODY const &body){  
  forall<typename POL::pol_l>(0, end_l, [&](int l){  
    forall<typename POL::pol_k>(0, end_k, [&](int k){
      forall<typename POL::pol_i>(0, end_i, [&](int i){
        forall<typename POL::pol_j>(0, end_j, [&](int j){        
          body(i,j,k,l);  
        });
      });
    });
  });
}

template<typename POL, typename BODY>
inline void forall4(LAYOUT_LKJI const, int end_i, int end_j, int end_k, int end_l, BODY const &body){  
  forall<typename POL::pol_l>(0, end_l, [&](int l){  
    forall<typename POL::pol_k>(0, end_k, [&](int k){      
      forall<typename POL::pol_j>(0, end_j, [&](int j){        
        forall<typename POL::pol_i>(0, end_i, [&](int i){
          body(i,j,k,l);  
        });
      });
    });
  });
}


template<typename POL, typename BODY>
inline void forall4(int end_i, int end_j, int end_k, int end_l, BODY const &body){
  typedef typename POL::layout L;
  forall4<POL, BODY>(L(), end_i, end_j, end_k, end_l, body);
}




#endif

