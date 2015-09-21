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

#include<Kripke/Kernel/Kernel_3d_ZGD.h>
#include<Kripke/Grid.h>
#include<Kripke/SubTVec.h>
#include<Domain/Layout.h>

Kernel_3d_ZGD::Kernel_3d_ZGD() :
  Kernel(NEST_ZGD)
{}

Kernel_3d_ZGD::~Kernel_3d_ZGD()
{}

Nesting_Order Kernel_3d_ZGD::nestingPsi(void) const {
  return NEST_ZGD;
}

Nesting_Order Kernel_3d_ZGD::nestingPhi(void) const {
  return NEST_ZGD;
}

Nesting_Order Kernel_3d_ZGD::nestingSigt(void) const {
  return NEST_DZG;
}

Nesting_Order Kernel_3d_ZGD::nestingEll(void) const {
  return NEST_ZGD;
}

Nesting_Order Kernel_3d_ZGD::nestingEllPlus(void) const {
  return NEST_ZDG;
}

Nesting_Order Kernel_3d_ZGD::nestingSigs(void) const {
  return NEST_ZGD;
}


/**
 * Add an isotropic source, with flux of 1, to every zone with Region 1
 * (or material 0).
 *
 * Since it's isotropic, we're just adding this to nm=0.
 */
void Kernel_3d_ZGD::source(Grid_Data *grid_data){
  // Loop over zoneset subdomains
  for(int zs = 0;zs < grid_data->num_zone_sets;++ zs){
  
    // get material mix information
    int sdom_id = grid_data->zs_to_sdomid[zs];
    Subdomain &sdom = grid_data->subdomains[sdom_id];

    // grab dimensions
    int num_mixed = sdom.mixed_to_zones.size();
    int num_zones = sdom.num_zones;
    int num_groups = grid_data->phi_out[zs]->groups;
    int num_moments = grid_data->total_num_moments;
    
    View3d<double, LAYOUT_KJI> phi_out(sdom.phi_out->ptr(), num_moments, num_groups, num_zones);
    View1d<int,    LAYOUT_I> const mixed_to_zones(&sdom.mixed_to_zones[0], 1);
    View1d<int,    LAYOUT_I> const mixed_material(&sdom.mixed_material[0], 1);
    View1d<double, LAYOUT_I> const mixed_fraction(&sdom.mixed_fraction[0], 1);

#ifdef KRIPKE_USE_OPENMP
#pragma omp parallel for
#endif
    for(int mix = 0;mix < num_mixed;++ mix){
      int zone = mixed_to_zones(mix);
      int material = mixed_material(mix);
      double fraction = mixed_fraction(mix);
        
      if(material == 0){        
        for(int g = 0;g < num_groups;++ g){
          phi_out(0, g, zone) += 1.0 * fraction;        
        }
      }
    }
  }
}


void Kernel_3d_ZGD::sweep(Subdomain *sdom) {
  int num_directions = sdom->num_directions;
  int num_groups = sdom->num_groups;
  int num_zones = sdom->num_zones;

  Directions *direction = sdom->directions;

  int local_imax = sdom->nzones[0];
  int local_jmax = sdom->nzones[1];
  int local_kmax = sdom->nzones[2];

  View1d<double, LAYOUT_I> const dx(&sdom->deltas[0][0], local_imax+2);
  View1d<double, LAYOUT_I> const dy(&sdom->deltas[1][0], local_jmax+2);
  View1d<double, LAYOUT_I> const dz(&sdom->deltas[2][0], local_kmax+2);
  
  View3d<double, LAYOUT_KJI> const rhs(sdom->rhs->ptr(), num_directions, num_groups, num_zones);
  View3d<double, LAYOUT_KJI> psi(sdom->psi->ptr(), num_directions, num_groups, num_zones);
  View2d<double, LAYOUT_JI>  const sigt(sdom->sigt->ptr(), num_groups, num_zones);

  int num_z_i = local_jmax * local_kmax;
  int num_z_j = local_imax * local_kmax;
  int num_z_k = local_imax * local_jmax;  
  
  View3d<double, LAYOUT_KJI> psi_lf(sdom->plane_data[0]->ptr(), num_directions, num_groups, num_z_i);
  View3d<double, LAYOUT_KJI> psi_fr(sdom->plane_data[1]->ptr(), num_directions, num_groups, num_z_j);
  View3d<double, LAYOUT_KJI> psi_bo(sdom->plane_data[2]->ptr(), num_directions, num_groups, num_z_k);
  
  Layout3d<LAYOUT_KJI> zone_layout(local_imax, local_jmax, local_kmax);
  Layout2d<LAYOUT_JI> i_layout(local_jmax, local_kmax);
  Layout2d<LAYOUT_JI> j_layout(local_imax, local_kmax);
  Layout2d<LAYOUT_JI> k_layout(local_imax, local_jmax);
  

  // All directions have same id,jd,kd, since these are all one Direction Set
  // So pull that information out now
  Grid_Sweep_Block const &extent = sdom->sweep_block;

  //  Perform transport sweep of the grid 1 cell at a time.
  for (int k = extent.start_k; k != extent.end_k; k += extent.inc_k) {
    for (int j = extent.start_j; j != extent.end_j; j += extent.inc_j) {
      for (int i = extent.start_i; i != extent.end_i; i += extent.inc_i) {
#ifdef KRIPKE_USE_OPENMP
#pragma omp parallel for
#endif
        for (int g = 0; g < num_groups; ++g) {
          for (int d = 0; d < num_directions; ++d) {
            double const xcos_dxi = 2.0 * direction[d].xcos / dx(i + 1);
            double const ycos_dyj = 2.0 * direction[d].ycos / dy(j + 1);
            double const zcos_dzk = 2.0 * direction[d].zcos / dz(k + 1);
            
            int const z_idx = zone_layout(i,j,k);            
            int const lf_idx = i_layout(j,k);
            int const fr_idx = j_layout(i,k);
            int const bo_idx = k_layout(i,j);

            /* Calculate new zonal flux */
            double const psi_d_g_z = (
                  rhs(d,g,z_idx)                  
                + psi_lf(d,g,lf_idx) * xcos_dxi
                + psi_fr(d,g,fr_idx) * ycos_dyj
                + psi_bo(d,g,bo_idx) * zcos_dzk)
                / (xcos_dxi + ycos_dyj + zcos_dzk + sigt(g,z_idx) );

            psi(d,g,z_idx) = psi_d_g_z;
            
            /* Apply diamond-difference relationships */
            psi_lf(d,g,lf_idx) = 2.0 * psi_d_g_z - psi_lf(d,g,lf_idx);
            psi_fr(d,g,fr_idx) = 2.0 * psi_d_g_z - psi_fr(d,g,fr_idx);
            psi_bo(d,g,bo_idx) = 2.0 * psi_d_g_z - psi_bo(d,g,bo_idx);
          }
        }
      }
    }
  }
}

