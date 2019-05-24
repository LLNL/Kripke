/*
 * NOTICE
 *
 * This work was produced at the Lawrence Livermore National Laboratory (LLNL)
 * under contract no. DE-AC-52-07NA27344 (Contract 44) between the U.S.
 * Department of Energy (DOE) and Lawrence Livermore National Security, LLC
 * (LLNS) for the operation of LLNL. The rights of the Federal Government are
 * reserved under Contract 44.
 *
 * DISCLAIMER
 *
 * This work was prepared as an account of work sponsored by an agency of the
 * United States Government. Neither the United States Government nor Lawrence
 * Livermore National Security, LLC nor any of their employees, makes any
 * warranty, express or implied, or assumes any liability or responsibility
 * for the accuracy, completeness, or usefulness of any information, apparatus,
 * product, or process disclosed, or represents that its use would not infringe
 * privately-owned rights. Reference herein to any specific commercial products,
 * process, or service by trade name, trademark, manufacturer or otherwise does
 * not necessarily constitute or imply its endorsement, recommendation, or
 * favoring by the United States Government or Lawrence Livermore National
 * Security, LLC. The views and opinions of authors expressed herein do not
 * necessarily state or reflect those of the United States Government or
 * Lawrence Livermore National Security, LLC, and shall not be used for
 * advertising or product endorsement purposes.
 *
 * NOTIFICATION OF COMMERCIAL USE
 *
 * Commercialization of this product is prohibited without notifying the
 * Department of Energy (DOE) or Lawrence Livermore National Security.
 */

#include <Kripke/Kernel.h>

#include <Kripke.h>
#include <Kripke/Arch/SweepSubdomains.h>
#include <Kripke/Arch/Population.h>
#include <Kripke/Core/Comm.h>
#include <Kripke/Timing.h>
#include <Kripke/VarTypes.h>

#include <Kripke/KSN_AnalyticalSolutions.hh>

using namespace Kripke;
using namespace Kripke::Core;

using namespace KSN_AnalyticalSolutions;

struct ErrorNormsSdom {

  template<typename AL>
  void operator()(AL al, 
                  Kripke::SdomId sdom_id,
		  Kripke::Core::DataStore &data_store,
                  double                  *max_nrm_ptr,
                  double                  *l1_nrm_ptr,
                  double                  *l2_nrm_ptr,
                  double                  *max_sol_nrm_ptr,
                  double                  *l1_sol_nrm_ptr,
                  double                  *l2_sol_nrm_ptr,
		  double                  *vol_grp) const
  {

    using PopulationPolicy = Kripke::Arch::Policy_Population<AL>; // use population's policy to only to grab the reduction policy, avoids needing new policy for errornorms
    using ReducePolicy = typename PopulationPolicy::ReducePolicy; 

    using ExecPolicy = typename Kripke::Arch::Policy_SweepSubdomains<AL>::ExecPolicy; // we can re-use this policy because this uses a similar loop structure

    auto sdom_al = getSdomAL(al, sdom_id);

    int num_directions = data_store.getVariable<Set>("Set/Direction").size(sdom_id);
    int num_groups = data_store.getVariable<Set>("Set/Group").size(sdom_id);
    int local_imax = data_store.getVariable<Set>("Set/ZoneI").size(sdom_id);
    int local_jmax = data_store.getVariable<Set>("Set/ZoneJ").size(sdom_id);
    int local_kmax = data_store.getVariable<Set>("Set/ZoneK").size(sdom_id);

    size_t global_i0 = data_store.getVariable<Set>("Set/ZoneI").lower(sdom_id);
    size_t global_j0 = data_store.getVariable<Set>("Set/ZoneJ").lower(sdom_id);
    size_t global_k0 = data_store.getVariable<Set>("Set/ZoneK").lower(sdom_id);

    auto xcos = sdom_al.getView(data_store.getVariable<Field_Direction2Double>("quadrature/xcos"));
    auto ycos = sdom_al.getView(data_store.getVariable<Field_Direction2Double>("quadrature/ycos"));
    auto zcos = sdom_al.getView(data_store.getVariable<Field_Direction2Double>("quadrature/zcos"));
    auto view_id = sdom_al.getView(data_store.getVariable<Field_Direction2Int>("quadrature/id"));
    auto view_jd = sdom_al.getView(data_store.getVariable<Field_Direction2Int>("quadrature/jd"));
    auto view_kd = sdom_al.getView(data_store.getVariable<Field_Direction2Int>("quadrature/kd"));

    auto dx = sdom_al.getView(data_store.getVariable<Field_ZoneI2Double>("dx"));
    auto dy = sdom_al.getView(data_store.getVariable<Field_ZoneJ2Double>("dy"));
    auto dz = sdom_al.getView(data_store.getVariable<Field_ZoneK2Double>("dz"));
    
    auto &field_w =         data_store.getVariable<Field_Direction2Double>("quadrature/w");
    auto &field_volume =    data_store.getVariable<Field_Zone2Double>("volume");
    auto w      = sdom_al.getView(field_w);
    auto volume = sdom_al.getView(field_volume);

    auto psi = sdom_al.getView(data_store.getVariable<Kripke::Field_Flux>("psi"));

    auto zone_layout = data_store.getVariable<ProductSet<3>>("Set/Zone").getLayout(sdom_id);

    ZoneI start_i(0);
    ZoneJ start_j(0);
    ZoneK start_k(0);

    ZoneI end_i(local_imax);
    ZoneJ end_j(local_jmax);
    ZoneK end_k(local_kmax);

    RAJA::ReduceMax<ReducePolicy, double> max_nrm_red(0.0);
    RAJA::ReduceSum<ReducePolicy, double> l1_nrm_red(0.0);
    RAJA::ReduceSum<ReducePolicy, double> l2_nrm_red(0.0);
    RAJA::ReduceMax<ReducePolicy, double> max_sol_nrm_red(0.0);
    RAJA::ReduceSum<ReducePolicy, double> l1_sol_nrm_red(0.0);
    RAJA::ReduceSum<ReducePolicy, double> l2_sol_nrm_red(0.0);
    RAJA::ReduceSum<ReducePolicy, double> vol_grp_red(0.0);

    Problem_3<double> ksn_problem;
    Direction d0{0};
    //int id = view_id(d0);
    //int jd = view_jd(d0);
    //int kd = view_kd(d0);
    
    RAJA::kernel<ExecPolicy>(
			     camp::make_tuple(
					      RAJA::TypedRangeSegment<Direction>(0, num_directions),
					      RAJA::TypedRangeSegment<Group>(0, num_groups),
					      RAJA::TypedRangeStrideSegment<ZoneK>(*start_k, *end_k, 1),
					      RAJA::TypedRangeStrideSegment<ZoneJ>(*start_j, *end_j, 1),
					      RAJA::TypedRangeStrideSegment<ZoneI>(*start_i, *end_i, 1)
					      

					      ),
			     KRIPKE_LAMBDA (Direction d, Group g, ZoneK k, ZoneJ j, ZoneI i) {


			       double xz = -60  + dx(i)*global_i0 + dx(i)*(double(*i)+0.5);
			       double yz = -100 + dy(j)*global_j0 + dy(j)*(double(*j)+0.5);
			       double zz = -60  + dz(k)*global_k0 + dz(k)*(double(*k)+0.5);

			       Zone z(zone_layout(*i, *j, *k));
			       double psi_exact = compute_exact_solution<double>(xz,yz,zz,
										 view_id(d)*xcos(d),view_jd(d)*ycos(d),view_kd(d)*zcos(d),
										 ksn_problem);
			       double err = std::abs(psi(d,g,z) - psi_exact);

						     
			       max_nrm_red.max(err);
			       l1_nrm_red += w(d) * volume(z) * err;
			       l2_nrm_red += w(d) * volume(z) * err * err;

			       max_sol_nrm_red.max(std::abs(psi(d,g,z)));
			       l1_sol_nrm_red += w(d) * volume(z) * std::abs(psi(d,g,z));
			       l2_sol_nrm_red += w(d) * volume(z) * std::abs(psi(d,g,z)) * std::abs(psi(d,g,z));

			       vol_grp_red += w(d) * volume(z);
			     }
			     );
    
    *max_nrm_ptr = fmax(*max_nrm_ptr,(double) max_nrm_red);
    *l1_nrm_ptr += (double) l1_nrm_red;
    *l2_nrm_ptr += (double) l2_nrm_red;

    *max_sol_nrm_ptr = fmax(*max_sol_nrm_ptr,(double) max_sol_nrm_red);
    *l1_sol_nrm_ptr += (double) l1_sol_nrm_red;
    *l2_sol_nrm_ptr += (double) l2_sol_nrm_red;

    *vol_grp += (double) vol_grp_red;
  }

};


/**
 * Returns the integral of Psi over all phase-space, to look at convergence
 */
void Kripke::Kernel::error_norms(Kripke::Core::DataStore &data_store)
{
  KRIPKE_TIMER(data_store, ErrorNorms);

  ArchLayoutV al_v = data_store.getVariable<ArchLayout>("al").al_v; 

  // sum up particles for psi and rhs
  double max_norm = 0.0;
  double l1_norm = 0.0;
  double l2_norm = 0.0;
  double max_sol_norm = 0.0;
  double l1_sol_norm = 0.0;
  double l2_sol_norm = 0.0;
  double vol_grp = 0.0;
  auto &field_psi = data_store.getVariable<Kripke::Field_Flux>("psi");
    
  for (Kripke::SdomId sdom_id : field_psi.getWorkList()){

    Kripke::dispatch(al_v, ErrorNormsSdom{}, sdom_id,
		     data_store,
                     &max_norm, &l1_norm, &l2_norm, &max_sol_norm, &l1_sol_norm, &l2_sol_norm, &vol_grp);
  }

  // reduce
  auto const &comm = data_store.getVariable<Kripke::Core::Comm>("comm");

  double vol = comm.allReduceSumDouble(vol_grp);
  double g_max_norm = comm.allReduceMaxDouble(max_norm);
  double g_max_sol_norm = comm.allReduceMaxDouble(max_sol_norm);

  double g_l1_norm = comm.allReduceSumDouble(l1_norm)/vol;
  double g_l2_norm = std::sqrt(comm.allReduceSumDouble(l2_norm)/vol);

  double g_l1_sol_norm = comm.allReduceSumDouble(l1_sol_norm)/vol;
  double g_l2_sol_norm = std::sqrt(comm.allReduceSumDouble(l2_sol_norm)/vol);

  if(comm.rank() == 0){
    printf("  errors(max, l2, l1) : %18.12e, %18.12e, %18.12e\n", g_max_norm, g_l2_norm, g_l1_norm);
    printf("  solution(max, l2, l1) : %18.12e, %18.12e, %18.12e\n", g_max_sol_norm, g_l2_sol_norm, g_l1_sol_norm);
    fflush(stdout);
  }
  
  return ;
}
