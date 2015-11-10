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

#ifdef _OPENMP

typedef RAJA::simd_exec seq_pol;
typedef RAJA::omp_parallel_for_exec omp_pol;
typedef RAJA::omp_for_nowait_exec omp_nowait;
typedef RAJA::omp_parallel_seq_exec omp_parallel_seq;

typedef RAJA::IndexSet::ExecPolicy<RAJA::seq_segit, RAJA::simd_exec> sweep_seq_pol;
typedef RAJA::IndexSet::ExecPolicy<RAJA::seq_segit, RAJA::omp_parallel_for_exec> sweep_omp_pol;
typedef RAJA::IndexSet::ExecPolicy<RAJA::seq_segit, RAJA::omp_parallel_seq_exec> sweep_parallel_seq_pol;


#else

typedef RAJA::simd_exec seq_pol;
typedef RAJA::simd_exec omp_pol;
typedef RAJA::simd_exec omp_nowait;
typedef RAJA::simd_exec omp_parallel_seq;

typedef RAJA::IndexSet::ExecPolicy<RAJA::seq_segit, RAJA::simd_exec> sweep_seq_pol;
typedef RAJA::IndexSet::ExecPolicy<RAJA::seq_segit, RAJA::simd_exec> sweep_omp_pol;
typedef RAJA::IndexSet::ExecPolicy<RAJA::seq_segit, RAJA::simd_exec> sweep_parallel_seq_pol;


#endif


// Include nested forall's
#include<Domain/Forall2.h>
#include<Domain/Forall3.h>
#include<Domain/Forall4.h>



// Subdomain loops
template<typename SubdomainPolicy, typename BODY>
RAJA_INLINE void forallSubdomains(Grid_Data *grid_data, BODY const &body){

  RAJA::forall<SubdomainPolicy>(
    RAJA::RangeSegment(0, grid_data->subdomains.size()),
    [=](int sdom_id){
      // get subdomain object
      Subdomain &sdom = grid_data->subdomains[sdom_id];

      body(sdom_id, sdom);
    });

}

// Loop over zoneset subdomains
template<typename SubdomainPolicy, typename BODY>
RAJA_INLINE void forallZoneSets(Grid_Data *grid_data, BODY const &body){

  RAJA::forall<SubdomainPolicy>(
    RAJA::RangeSegment(0, grid_data->num_zone_sets),
    [=](int zs){
      // get material mix information
      int sdom_id = grid_data->zs_to_sdomid[zs];
      Subdomain &sdom = grid_data->subdomains[sdom_id];

      body(zs, sdom_id, sdom);
    });

}





// Foreward decl
template<typename TAG, typename BODY>
RAJA_INLINE void executeScope(BODY const &body);

// Foreward decl
template<typename TAG, typename T1, typename L1, typename... REST>
RAJA_INLINE void executeScope(T1 t1, L1 const &l1, REST... rest);


template<typename T, typename TAG>
struct ExecuteIfMatches{
    template<typename BODY>
    RAJA_INLINE void match(BODY const &body) const {
      // nop
    }

    template<typename... REST>
    RAJA_INLINE void not_match(REST... rest) const {
      executeScope<TAG, REST...>(rest...);
    }
};


template<typename T>
struct ExecuteIfMatches<T, T>{
    template<typename BODY>
    RAJA_INLINE void match(BODY const &body) const {
      body(T());
    }

    template<typename... REST>
    RAJA_INLINE void not_match(REST... rest) const {
      // NOP
    }
};


/**
 * Termination case: labmda is always executed
 *
 * This handles the default policyScope lambd
 */
template<typename TAG, typename BODY>
RAJA_INLINE void executeScope(BODY const &body){
  body(TAG());
}

/**
 * Executes the L1 lambda with parameter type T1 if the types T1 and TAG are the same.
 *
 * If they are not the same, it recursively calls itself with the remaining types and arguments in REST... and rest...
 */
template<typename TAG, typename T1, typename L1, typename... REST>
RAJA_INLINE void executeScope(T1 t1, L1 const &l1, REST... rest){
  ExecuteIfMatches<T1, TAG> test;
  test.match(l1); // execute l1(T1()) if it matches
  test.not_match(rest...); // otherwise, try the remaining specializations
}


/**
 * Converts run-time nesting to static type for a given lambda scope.
 * This requires C++14 for polymorphic lambdas.
 *
 * First argument is the nesting order variable.
 * Followed by pairs of specializations [nest_tag(), labmda(nest_tag)],
 * Terminated by a single default lambda(auto nest_tag).
 *
 * The specializations are executed if any of their types match the nesting order tag.
 * The default is executed if none of the specializations match.
 *
 */
template<typename... REST>
RAJA_INLINE void policyScope(Nesting_Order nest, REST... rest){
  switch(nest){
    case NEST_DGZ: executeScope<NEST_DGZ_T>(rest...); break;
    case NEST_DZG: executeScope<NEST_DZG_T>(rest...); break;
    case NEST_GDZ: executeScope<NEST_GDZ_T>(rest...); break;
    case NEST_GZD: executeScope<NEST_GZD_T>(rest...); break;
    case NEST_ZDG: executeScope<NEST_ZDG_T>(rest...); break;
    case NEST_ZGD: executeScope<NEST_ZGD_T>(rest...); break;
  }
}


#endif
