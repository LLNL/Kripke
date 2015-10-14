#ifndef DOMAIN_FORALL_H__
#define DOMAIN_FORALL_H__

#include<Kripke/Subdomain.h>
#include<Domain/Layout.h>
#include<Domain/View.h>

#include<RAJA/RAJA.hxx>

//#define RAJA_INLINE __attribute__((always_inline))


//#define RAJA_LAMBDA [&]
#define RAJA_LAMBDA [=]
//#define RAJA_LAMBDA [=] __device__

typedef RAJA::simd_exec seq_pol;
typedef RAJA::omp_parallel_for_exec omp_pol;
typedef RAJA::omp_for_nowait_exec omp_nowait;

typedef RAJA::IndexSet::ExecPolicy<RAJA::seq_segit, RAJA::simd_exec> sweep_seq_pol;
typedef RAJA::IndexSet::ExecPolicy<RAJA::seq_segit, RAJA::omp_parallel_for_exec> sweep_omp_pol;
typedef RAJA::IndexSet::ExecPolicy<RAJA::seq_segit, RAJA::omp_parallel_seq_exec> sweep_parallel_seq_pol;

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
 * This requires C++14 for polymorphic lambdas
 */
template<typename BODY>
inline void policyScope(Nesting_Order nesting_order, BODY const &body){
  switch(nesting_order){
    case NEST_DGZ: body(NEST_DGZ_T()); break;
    case NEST_DZG: body(NEST_DZG_T()); break;
    case NEST_GDZ: body(NEST_GDZ_T()); break;
    case NEST_GZD: body(NEST_GZD_T()); break;
    case NEST_ZDG: body(NEST_ZDG_T()); break;
    case NEST_ZGD: body(NEST_ZGD_T()); break;
  }
}



#endif
