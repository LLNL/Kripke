#ifndef DOMAIN_FORALL_H__
#define DOMAIN_FORALL_H__

#include<Kripke/Subdomain.h>
#include<Domain/Layout.h>
#include<Domain/View.h>

#include<RAJA/RAJA.hxx>

//#define RAJA_INLINE __attribute__((always_inline))


#define RAJA_LAMBDA [&]
//#define RAJA_LAMBDA [=]

typedef RAJA::simd_exec seq_pol;
//typedef RAJA::IndexSet::ExecPolicy<RAJA::seq_segit, RAJA::seq_exec> seq_pol;
//typedef RAJA::IndexSet::ExecPolicy<RAJA::seq_exec, RAJA::omp_parallel_segit> omp_pol;

typedef seq_pol omp_pol;

typedef RAJA::IndexSet::ExecPolicy<RAJA::seq_segit, RAJA::simd_exec> sweep_seq_pol;


// Include nested forall's
#include<Domain/Forall2.h>
#include<Domain/Forall3.h>
#include<Domain/Forall4.h>


// Subdomain loops
template<typename BODY>
RAJA_INLINE void forallSubdomains(Grid_Data *grid_data, BODY const &body){

  int num_subdomains = grid_data->subdomains.size();
  for (int sdom_id = 0; sdom_id < num_subdomains; ++ sdom_id){
    Subdomain &sdom = grid_data->subdomains[sdom_id];

    body(sdom_id, sdom);
  }

}

// Loop over zoneset subdomains
template<typename BODY>
RAJA_INLINE void forallZoneSets(Grid_Data *grid_data, BODY const &body){

  for(int zs = 0;zs < grid_data->num_zone_sets;++ zs){
    // get material mix information
    int sdom_id = grid_data->zs_to_sdomid[zs];
    Subdomain &sdom = grid_data->subdomains[sdom_id];

    body(zs, sdom_id, sdom);
  }

}




/**
 * Converts run-time nesting to static type for a given lambda scope.
 * This could be done w/o a macro in C++14 with polymorphic lambdas
 */
//#if 1 // set to 0 to use polymo

template<typename BODY>
inline void policyScopeT(Nesting_Order nesting_order, BODY const &body){
  switch(nesting_order){
    case NEST_DGZ: body(NEST_DGZ_T()); break;
    case NEST_DZG: body(NEST_DZG_T()); break;
    case NEST_GDZ: body(NEST_GDZ_T()); break;
    case NEST_GZD: body(NEST_GZD_T()); break;
    case NEST_ZDG: body(NEST_ZDG_T()); break;
    case NEST_ZGD: body(NEST_ZGD_T()); break;
  }
}


//#else

#define policyScope(GNT_NEST__, GNT_TYPE__, GNT_LAMBDA__) \
{\
  switch(GNT_NEST__){\
    case NEST_DGZ: {\
      typedef NEST_DGZ_T GNT_TYPE__;\
      GNT_LAMBDA__();\
    } break; \
    case NEST_DZG: {\
      typedef NEST_DZG_T GNT_TYPE__;\
      GNT_LAMBDA__();\
    } break; \
    case NEST_GDZ: {\
      typedef NEST_GDZ_T GNT_TYPE__;\
      GNT_LAMBDA__();\
    } break; \
    case NEST_GZD: {\
      typedef NEST_GZD_T GNT_TYPE__;\
      GNT_LAMBDA__();\
    } break; \
    case NEST_ZDG: {\
      typedef NEST_ZDG_T GNT_TYPE__;\
      GNT_LAMBDA__();\
    } break; \
    case NEST_ZGD: {\
      typedef NEST_ZGD_T GNT_TYPE__;\
      GNT_LAMBDA__();\
    } break; \
  }\
}

//#endif

#if 1

#define BEGIN_POLICY_SCOPE(NEST, TYPE) policyScopeT(NEST, [&](auto POLTYPE){\
                                         typedef decltype(POLTYPE) TYPE;
#define   END_POLICY_SCOPE()             });

#else

#define BEGIN_POLICY_SCOPE(NEST, TYPE) policyScope(NEST, TYPE, [&](){
#define   END_POLICY_SCOPE()             });


#endif


#endif

