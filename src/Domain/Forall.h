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

struct seq2_pol {
  typedef seq_pol pol_i;
  typedef seq_pol pol_j;
};
/*
template<typename POL, typename BODY>
inline void forall2(int start_i, int end_i, int start_j, int end_j, BODY const &body){
  forall<typename POL::pol_i>(start_i, end_i, [&](int i){
    forall<typename POL::pol_j>(start_j, end_j, [&](int j){
      body(i,j);
    });
  });
}*/


template<typename POL, typename BODY>
inline void forall4(LAYOUT_IJKL_t const, int end_i, int end_j, int end_k, int end_l, BODY const &body){
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
inline void forall4(LAYOUT_IJLK_t const, int end_i, int end_j, int end_k, int end_l, BODY const &body){
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
inline void forall4(int end_i, int end_j, int end_k, int end_l, BODY const &body){
  typedef typename POL::layout_t L;
  forall4<POL, BODY>(L(), end_i, end_j, end_k, end_l, body);
}

/*
template<typename POL, typename BODY>
inline void forall4(int end_i, int end_j, int end_k, int end_l, BODY const &body){
  switch(POL::layout){
    case LAYOUT_IJKL:
      forall<typename POL::pol_i>(0, end_i, [&](int i){
        forall<typename POL::pol_j>(0, end_j, [&](int j){
          forall<typename POL::pol_k>(0, end_k, [&](int k){
            forall<typename POL::pol_l>(0, end_l, [&](int l){
              body(i,j,k,l);  
            });
          });
        });
      });
      break;
      
    case LAYOUT_IJLK:
      forall<typename POL::pol_i>(0, end_i, [&](int i){
        forall<typename POL::pol_j>(0, end_j, [&](int j){          
          forall<typename POL::pol_l>(0, end_l, [&](int l){
            forall<typename POL::pol_k>(0, end_k, [&](int k){
              body(i,j,k,l);  
            });
          });
        });
      });
      break;
  }
}*/

#endif

