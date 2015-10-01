#ifndef DOMAIN_FORALL_H__
#define DOMAIN_FORALL_H__

#include<Kripke/Subdomain.h>
#include<Domain/Layout.h>
#include<Domain/View.h>

#include<RAJA/RAJA.hxx>

//#define RAJA_INLINE __attribute__((always_inline))

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


#endif

