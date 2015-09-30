#ifndef DOMAIN_FORALL_H__
#define DOMAIN_FORALL_H__

#include<Domain/Layout.h>
#include<Domain/View.h>

struct seq_pol{};
struct omp_pol{};


template<typename POL, typename BODY>
inline void forall(int start, int end, BODY const &body){
  forall(POL(), start, end, body);
}

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


// Include nested forall's
#include<Domain/Forall2.h>
#include<Domain/Forall4.h>


#endif

