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
#include <Kripke/Timing.h>
#include <Kripke/VarTypes.h>

using namespace Kripke;
using namespace Kripke::Core;



struct SweepSdom {

  template<typename AL>
  RAJA_INLINE
  void operator()(AL, Kripke::Core::DataStore &data_store,
                  Kripke::SdomId              sdom_id) const
  {

  
    int num_directions = data_store.getVariable<Set>("Set/Direction").size(sdom_id);
    int num_groups = data_store.getVariable<Set>("Set/Group").size(sdom_id);
    int local_imax = data_store.getVariable<Set>("Set/ZoneI").size(sdom_id);
    int local_jmax = data_store.getVariable<Set>("Set/ZoneJ").size(sdom_id);
    int local_kmax = data_store.getVariable<Set>("Set/ZoneK").size(sdom_id);

    auto xcos = data_store.getVariable<Field_Direction2Double>("quadrature/xcos").getViewAL<AL>(sdom_id);
    auto ycos = data_store.getVariable<Field_Direction2Double>("quadrature/ycos").getViewAL<AL>(sdom_id);
    auto zcos = data_store.getVariable<Field_Direction2Double>("quadrature/zcos").getViewAL<AL>(sdom_id);
    auto view_id = data_store.getVariable<Field_Direction2Int>("quadrature/id").getViewAL<AL>(sdom_id);
    auto view_jd = data_store.getVariable<Field_Direction2Int>("quadrature/jd").getViewAL<AL>(sdom_id);
    auto view_kd = data_store.getVariable<Field_Direction2Int>("quadrature/kd").getViewAL<AL>(sdom_id);

    auto dx = data_store.getVariable<Field_ZoneI2Double>("dx").getViewAL<AL>(sdom_id);
    auto dy = data_store.getVariable<Field_ZoneJ2Double>("dy").getViewAL<AL>(sdom_id);
    auto dz = data_store.getVariable<Field_ZoneK2Double>("dz").getViewAL<AL>(sdom_id);

    auto sigt = data_store.getVariable<Kripke::Field_SigmaTZonal>("sigt_zonal").getViewAL<AL>(sdom_id);
    auto psi = data_store.getVariable<Kripke::Field_Flux>("psi").getViewAL<AL>(sdom_id);
    auto rhs = data_store.getVariable<Kripke::Field_Flux>("rhs").getViewAL<AL>(sdom_id);

    auto psi_lf = data_store.getVariable<Field_IPlane>("i_plane").getViewAL<AL>(sdom_id);
    auto psi_fr = data_store.getVariable<Field_JPlane>("j_plane").getViewAL<AL>(sdom_id);
    auto psi_bo = data_store.getVariable<Field_KPlane>("k_plane").getViewAL<AL>(sdom_id);

    // Assumption: all directions in this sdom have same mesh traversal
    Direction d0{0};
    int id = view_id(d0);
    int jd = view_jd(d0);
    int kd = view_kd(d0);

    ZoneI start_i((id>0) ? 0 : local_imax-1);
    ZoneJ start_j((jd>0) ? 0 : local_jmax-1);
    ZoneK start_k((kd>0) ? 0 : local_kmax-1);

    ZoneI end_i((id>0) ? local_imax : -1);
    ZoneJ end_j((jd>0) ? local_jmax : -1);
    ZoneK end_k((kd>0) ? local_kmax : -1);

    auto zone_layout = data_store.getVariable<ProductSet<3>>("Set/Zone").getLayout(sdom_id);

    RAJA::kernel<Kripke::Arch::Policy_SweepSubdomains>(
        camp::make_tuple(
            RAJA::TypedRangeSegment<Direction>(0, num_directions),
            RAJA::TypedRangeSegment<Group>(0, num_groups),
            RAJA::TypedRangeStrideSegment<ZoneK>(*start_k, *end_k, kd),
            RAJA::TypedRangeStrideSegment<ZoneJ>(*start_j, *end_j, jd),
            RAJA::TypedRangeStrideSegment<ZoneI>(*start_i, *end_i, id)


    ),
        KRIPKE_LAMBDA (Direction d, Group g, ZoneK k, ZoneJ j, ZoneI i) {

            double xcos_dxi = 2.0 * xcos(d) / dx(i);
            double ycos_dyj = 2.0 * ycos(d) / dy(j);
            double zcos_dzk = 2.0 * zcos(d) / dz(k);

            Zone z(zone_layout(*i, *j, *k));

            /* Calculate new zonal flux */
            double psi_d_g_z = (rhs(d,g,z)
                + psi_lf(d, g, j, k) * xcos_dxi
                + psi_fr(d, g, i, k) * ycos_dyj
                + psi_bo(d, g, i, j) * zcos_dzk)
                / (xcos_dxi + ycos_dyj + zcos_dzk + sigt(g, z));

            psi(d, g, z) = psi_d_g_z;

            /* Apply diamond-difference relationships */
            psi_lf(d, g, j, k) = 2.0 * psi_d_g_z - psi_lf(d, g, j, k);
            psi_fr(d, g, i, k) = 2.0 * psi_d_g_z - psi_fr(d, g, i, k);
            psi_bo(d, g, i, j) = 2.0 * psi_d_g_z - psi_bo(d, g, i, j);

        }
    );
  }
};

void Kripke::Kernel::sweepSubdomain(Kripke::Core::DataStore &data_store,
    Kripke::SdomId sdom_id)
{
  KRIPKE_TIMER(data_store, SweepSubdomain);

  Kripke::Kernel::dispatch(data_store, sdom_id, SweepSdom{});

}
