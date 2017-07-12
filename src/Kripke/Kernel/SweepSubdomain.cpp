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
#include <Kripke/Grid.h>
#include <Kripke/SubTVec.h>
#include <Kripke/Timing.h>
#include <Kripke/VarTypes.h>


// Macros for offsets with fluxes on cell faces 
#define I_PLANE_INDEX(j, k) ((k)*(local_jmax) + (j))
#define J_PLANE_INDEX(i, k) ((k)*(local_imax) + (i))
#define K_PLANE_INDEX(i, j) ((j)*(local_imax) + (i))
#define Zonal_INDEX(i, j, k) ((i) + (local_imax)*(j) \
  + (local_imax)*(local_jmax)*(k))

void Kripke::Kernel::sweepSubdomain(Kripke::DataStore &data_store,
    Kripke::SdomId sdom_id)
{

  KRIPKE_TIMER(data_store, SweepSubdomain);

  auto &field_psi = data_store.getVariable<Kripke::Field_Flux>("psi");
  auto &field_rhs = data_store.getVariable<Kripke::Field_Flux>("rhs");

  Grid_Data *grid_data = &data_store.getVariable<Grid_Data>("grid_data");
  Subdomain *sdom = &(grid_data->subdomains[*sdom_id]);

  int num_directions = sdom->num_directions;
  int num_groups = sdom->num_groups;
  int num_zones = sdom->num_zones;

  Kripke::QuadraturePoint *direction = sdom->directions;

  int local_imax = sdom->nzones[0];
  int local_jmax = sdom->nzones[1];
  int local_kmax = sdom->nzones[2];

  double const * KRESTRICT dx = &sdom->deltas[0][0];
  double const * KRESTRICT dy = &sdom->deltas[1][0];
  double const * KRESTRICT dz = &sdom->deltas[2][0];
  
  double const * KRESTRICT sigt = sdom->sigt->ptr();
  double       * KRESTRICT psi = field_psi.getData(sdom_id);
  double const * KRESTRICT rhs = field_rhs.getData(sdom_id);

  double * KRESTRICT psi_lf = sdom->plane_data[0]->ptr();
  double * KRESTRICT psi_fr = sdom->plane_data[1]->ptr();
  double * KRESTRICT psi_bo = sdom->plane_data[2]->ptr();
  
  int num_gz = num_groups * num_zones;
  int num_gz_i = local_jmax * local_kmax * num_groups;
  int num_gz_j = local_imax * local_kmax * num_groups;
  int num_gz_k = local_imax * local_jmax * num_groups;
  int num_z_i = local_jmax * local_kmax;
  int num_z_j = local_imax * local_kmax;
  int num_z_k = local_imax * local_jmax;

  // All directions have same id,jd,kd, since these are all one Direction Set
  // So pull that information out now
  Grid_Sweep_Block const &extent = sdom->sweep_block;

  for (int d = 0; d < num_directions; ++d) {
    double xcos = 2.0 * direction[d].xcos;
    double ycos = 2.0 * direction[d].ycos;
    double zcos = 2.0 * direction[d].zcos;
    
    double       * KRESTRICT psi_d  = psi  + d*num_gz;
    double const * KRESTRICT rhs_d  = rhs  + d*num_gz;

    double       * KRESTRICT psi_lf_d = psi_lf + d*num_gz_i;
    double       * KRESTRICT psi_fr_d = psi_fr + d*num_gz_j;
    double       * KRESTRICT psi_bo_d = psi_bo + d*num_gz_k;

    for (int g = 0; g < num_groups; ++g) {
      double const * KRESTRICT sigt_g  = sigt + g*num_zones;
      double       * KRESTRICT psi_d_g = psi_d + g*num_zones;
      double const * KRESTRICT rhs_d_g = rhs_d + g*num_zones;
      
      double       * KRESTRICT psi_lf_d_g = psi_lf_d + g*num_z_i;
      double       * KRESTRICT psi_fr_d_g = psi_fr_d + g*num_z_j;
      double       * KRESTRICT psi_bo_d_g = psi_bo_d + g*num_z_k;

      for (int k = extent.start_k; k != extent.end_k; k += extent.inc_k) {       
        double zcos_dzk = zcos / dz[k + 1];
        
        for (int j = extent.start_j; j != extent.end_j; j += extent.inc_j) {
          double ycos_dyj = ycos / dy[j + 1];
          
          for (int i = extent.start_i; i != extent.end_i; i += extent.inc_i) {
            double xcos_dxi = xcos / dx[i + 1];
            
            int z_idx = Zonal_INDEX(i, j, k);
            int z_i = I_PLANE_INDEX(j, k);
            int z_j = J_PLANE_INDEX(i, k);
            int z_k = K_PLANE_INDEX(i, j);

            /* Calculate new zonal flux */
            double psi_d_g_z = (rhs_d_g[z_idx]
                + psi_lf_d_g[z_i] * xcos_dxi
                + psi_fr_d_g[z_j] * ycos_dyj
                + psi_bo_d_g[z_k] * zcos_dzk)
                / (xcos_dxi + ycos_dyj + zcos_dzk + sigt_g[z_idx]);

            psi_d_g[z_idx] = psi_d_g_z;
            
            /* Apply diamond-difference relationships */
            psi_lf_d_g[z_i] = 2.0 * psi_d_g_z - psi_lf_d_g[z_i];
            psi_fr_d_g[z_j] = 2.0 * psi_d_g_z - psi_fr_d_g[z_j];
            psi_bo_d_g[z_k] = 2.0 * psi_d_g_z - psi_bo_d_g[z_k];
          }
        }
      }
    } // group
  } // direction

}


