#ifndef DOMAIN_FORALL_H__
#define DOMAIN_FORALL_H__

#include<Domain/Layout.h>
#include<Domain/View.h>

#include<RAJA/RAJA.hxx>

typedef RAJA::seq_exec seq_pol;
//typedef RAJA::IndexSet::ExecPolicy<RAJA::seq_segit, RAJA::seq_exec> seq_pol;
//typedef RAJA::IndexSet::ExecPolicy<RAJA::seq_exec, RAJA::omp_parallel_segit> omp_pol;

typedef seq_pol omp_pol;


// Include nested forall's
#include<Domain/Forall2.h>
#include<Domain/Forall3.h>
#include<Domain/Forall4.h>

#endif

