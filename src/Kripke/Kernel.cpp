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

#include<Kripke/Kernel.h>
#include<Kripke/Grid.h>
#include<Kripke/SubTVec.h>

#include<Domain/View.h>
#include<Domain/Forall.h>

#include<RAJA/RAJA.hxx>

#include<Kripke/Kernel/Kernel_3d_GDZ.h>
#include<Kripke/Kernel/Kernel_3d_DGZ.h>
#include<Kripke/Kernel/Kernel_3d_ZDG.h>
#include<Kripke/Kernel/Kernel_3d_DZG.h>
#include<Kripke/Kernel/Kernel_3d_ZGD.h>
#include<Kripke/Kernel/Kernel_3d_GZD.h>

#include<Kripke/Kernel/VariablePolicy.h>

/**
 * Factory to create a kernel object for the specified nesting
 */
Kernel *createKernel(Nesting_Order nest, int num_dims){
  if(num_dims == 3){
    switch(nest){
    case NEST_GDZ:
      return new Kernel_3d_GDZ();
    case NEST_DGZ:
      return new Kernel_3d_DGZ();
    case NEST_ZDG:
      return new Kernel_3d_ZDG();
    case NEST_DZG:
      return new Kernel_3d_DZG();
    case NEST_ZGD:
      return new Kernel_3d_ZGD();
    case NEST_GZD:
      return new Kernel_3d_GZD();
    }
  }

  MPI_Abort(MPI_COMM_WORLD, 1);
  return NULL;
}


Kernel::Kernel(Nesting_Order nest) :
  nesting_order(nest)
{}

Kernel::~Kernel(){
}




#include<Kripke/Kernel/LTimesPolicy.h>
void Kernel::LTimes(Grid_Data *grid_data) { 
  // Outer parameters
  int num_moments = grid_data->total_num_moments;

  // Zero Phi
  forallZoneSets(grid_data, [=](int zs, int sdom_id, Subdomain &sdom){
    sdom.phi->clear(0.0);
  });

  // Loop over Subdomains
  forallSubdomains(grid_data, [=](int sdom_id, Subdomain &sdom){
  
    // Get dimensioning
    int num_zones = sdom.num_zones;
    int num_groups = sdom.phi->groups;
    int num_local_groups = sdom.num_groups;
    int group0 = sdom.group0;
    int num_local_directions = sdom.num_directions;

    policyScope(nesting_order, nest_type, [&](){
      typedef LayoutPolicy<nest_type> LAYOUT;
      typedef ViewPolicy<LAYOUT> VIEW;
                 
      // Get pointers    
      typename VIEW::Psi psi(sdom.psi->ptr(), num_local_directions, num_local_groups, num_zones);
      typename VIEW::Phi phi(sdom.phi->ptr(), num_moments, num_groups, num_zones);
      typename VIEW::Ell ell(sdom.ell->ptr(), num_local_directions, num_moments);
            
      forall4<LTimesPolicy<nest_type> >(
        RAJA::RangeSegment(0, num_moments),
        RAJA::RangeSegment(0, num_local_directions),
        RAJA::RangeSegment(0, num_local_groups),
        RAJA::RangeSegment(0, num_zones),
        [=](int nm, int d, int g, int z){
 
          phi(nm, g+group0, z) += ell(d,nm) * psi(d,g,z);
 
        }); // forall
    });// policy
  }); // sdom

}

#include<Kripke/Kernel/LPlusTimesPolicy.h>
void Kernel::LPlusTimes(Grid_Data *grid_data) {

  // Outer parameters
  int num_moments = grid_data->total_num_moments;

  // Loop over Subdomains
  forallSubdomains(grid_data, [=](int sdom_id, Subdomain &sdom){
    sdom.rhs->clear(0.0);
  });

  forallSubdomains(grid_data, [=](int sdom_id, Subdomain &sdom){

    // Get dimensioning
    int num_zones = sdom.num_zones;
    int num_local_groups = sdom.num_groups;
    int num_groups = sdom.phi_out->groups;
    int group0 = sdom.group0;
    int num_local_directions = sdom.num_directions;

    policyScope(nesting_order, nest_type, ([&](){
      typedef LayoutPolicy<nest_type> LAYOUT;
      typedef ViewPolicy<LAYOUT> VIEW;      
      // Get pointers
      VIEW::Psi     rhs(sdom.rhs->ptr(), num_local_directions, num_local_groups, num_zones);
      VIEW::Phi     phi_out(sdom.phi_out->ptr(), num_moments, num_groups, num_zones);
      VIEW::EllPlus ell_plus(sdom.ell_plus->ptr(), num_local_directions, num_moments);
      
      forall4<LPlusTimesPolicy<nest_type> >(
        RAJA::RangeSegment(0, num_moments),
        RAJA::RangeSegment(0, num_local_directions),
        RAJA::RangeSegment(0, num_local_groups),
        RAJA::RangeSegment(0, num_zones),
        [=](int nm, int d, int g, int z){
 
          rhs(d,g,z) += ell_plus(d,nm) * phi_out(nm,g+group0,z);
 
        }); // forall
    })); // policy
  }); // sdom
  
}


/**
  Compute scattering source term phi_out from flux moments in phi.
  phi_out(gp,z,nm) = sum_g { sigs(g, n, gp) * phi(g,z,nm) }
*/
#include<Kripke/Kernel/ScatteringPolicy.h>
void Kernel::scattering(Grid_Data *grid_data){

  // Zero out source terms
  forallZoneSets(grid_data, [=](int zs, int sdom_id, Subdomain &sdom){
    sdom.phi_out->clear(0.0);
  });
  
  // Loop over zoneset subdomains
  forallZoneSets(grid_data, [=](int zs, int sdom_id, Subdomain &sdom){

    // grab dimensions
    int num_zones = sdom.num_zones;
    int num_mixed = sdom.mixed_to_zones.size();
    int num_groups = grid_data->phi_out[zs]->groups;
    int num_moments = grid_data->total_num_moments;
    int legendre_order = grid_data->legendre_order;

    policyScope(nesting_order, nest_type, ([&](){
      typedef LayoutPolicy<nest_type> LAYOUT;
      typedef ViewPolicy<LAYOUT> VIEW;

      VIEW::Phi  phi_out(sdom.phi_out->ptr(), num_moments, num_groups, num_zones);
      VIEW::Phi  const phi(sdom.phi->ptr(), num_moments, num_groups, num_zones);  
      VIEW::SigS const sigs(grid_data->sigs->ptr(), legendre_order+1, num_groups, num_groups, 3);

      View1d<int, PERM_I>    const mixed_to_zones(&sdom.mixed_to_zones[0], 1);
      View1d<int, PERM_I>    const mixed_material(&sdom.mixed_material[0], 1);
      View1d<double, PERM_I> const mixed_fraction(&sdom.mixed_fraction[0], 1);
      View1d<int, PERM_I>    const moment_to_coeff(&grid_data->moment_to_coeff[0], 1);
          
      forall4<ScatteringPolicy<nest_type> >(
        RAJA::RangeSegment(0, num_moments),
        RAJA::RangeSegment(0, num_groups),
        RAJA::RangeSegment(0, num_groups),
        RAJA::RangeSegment(0, num_mixed),
        [=](int nm, int g, int gp, int mix){
        
          int n = moment_to_coeff(nm);
          int zone = mixed_to_zones(mix);
          int material = mixed_material(mix);
          double fraction = mixed_fraction(mix);

          phi_out(nm, gp, zone) += 
            sigs(n, g, gp, material) * phi(nm, g, zone) * fraction;                     
                             
        });  // forall
    })); // policy
  }); // zonesets
  
}

  
/**
 * Add an isotropic source, with flux of 1, to every zone with Region 1
 * (or material 0).
 *
 * Since it's isotropic, we're just adding this to nm=0.
 */
#include<Kripke/Kernel/SourcePolicy.h>
void Kernel::source(Grid_Data *grid_data){

  // Loop over zoneset subdomains
  forallZoneSets(grid_data, [=](int zs, int sdom_id, Subdomain &sdom){
       
    // grab dimensions
    int num_mixed = sdom.mixed_to_zones.size();
    int num_zones = sdom.num_zones;
    int num_groups = sdom.phi_out->groups;
    int num_moments = grid_data->total_num_moments;
    
    policyScope(nesting_order, nest_type, ([&](){
      typedef LayoutPolicy<nest_type> LAYOUT;
      typedef ViewPolicy<LAYOUT> VIEW;
    
      VIEW::Phi phi_out(sdom.phi_out->ptr(), num_moments, num_groups, num_zones);
      
      View1d<int,    PERM_I> const mixed_to_zones(&sdom.mixed_to_zones[0], 1);
      View1d<int,    PERM_I> const mixed_material(&sdom.mixed_material[0], 1);
      View1d<double, PERM_I> const mixed_fraction(&sdom.mixed_fraction[0], 1);

      forall2<SourcePolicy<nest_type> >(
        RAJA::RangeSegment(0, num_groups),
        RAJA::RangeSegment(0, num_mixed),
        [=](int g, int mix){
          int zone = mixed_to_zones(mix);
          int material = mixed_material(mix);
          double fraction = mixed_fraction(mix);

          if(material == 0){
            phi_out(0, g, zone) += 1.0 * fraction;
          }
        }); // forall
        
    })); // policy
  });// zonesets
}


#include<Kripke/Kernel/SweepPolicy.h>
void Kernel::sweep(Subdomain *sdom) {
  int num_directions = sdom->num_directions;
  int num_groups = sdom->num_groups;
  int num_zones = sdom->num_zones;

  Directions *direction = sdom->directions;

  int local_imax = sdom->nzones[0];
  int local_jmax = sdom->nzones[1];
  int local_kmax = sdom->nzones[2];
  
  int num_z_i = local_jmax * local_kmax;
  int num_z_j = local_imax * local_kmax;
  int num_z_k = local_imax * local_jmax;
    
  policyScope(nesting_order, nest_type, ([&](){
    typedef LayoutPolicy<nest_type> LAYOUT;
    typedef ViewPolicy<LAYOUT> VIEW;

    VIEW::Psi const rhs(sdom->rhs->ptr(), num_directions, num_groups, num_zones);
    VIEW::Psi psi(sdom->psi->ptr(), num_directions, num_groups, num_zones);
    VIEW::SigT const sigt(sdom->sigt->ptr(), num_groups, num_zones);
    
    View1d<double, PERM_I> const dx(&sdom->deltas[0][0], local_imax+2);
    View1d<double, PERM_I> const dy(&sdom->deltas[1][0], local_jmax+2);
    View1d<double, PERM_I> const dz(&sdom->deltas[2][0], local_kmax+2);
            
    VIEW::Face psi_lf(sdom->plane_data[0]->ptr(), num_directions, num_groups, num_z_i);
    VIEW::Face psi_fr(sdom->plane_data[1]->ptr(), num_directions, num_groups, num_z_j);
    VIEW::Face psi_bo(sdom->plane_data[2]->ptr(), num_directions, num_groups, num_z_k);
    
    LAYOUT::Zone zone_layout(local_imax, local_jmax, local_kmax);
    LAYOUT::Face i_layout(local_jmax, local_kmax);
    LAYOUT::Face j_layout(local_imax, local_kmax);
    LAYOUT::Face k_layout(local_imax, local_jmax);

    // All directions have same id,jd,kd, since these are all one Direction Set
    // So pull that information out now
    Grid_Sweep_Block const &extent = sdom->sweep_block;

    // TODO: We really need the RAJA::IndexSet concept to deal with zone iteration pattern
    forall3<SweepPolicy<nest_type> >(
      RAJA::RangeSegment(0, num_directions),
      RAJA::RangeSegment(0, num_groups),
      extent.indexset_sweep,
      [=](int d, int g, int zone_idx){

        int i = extent.idx_to_i[zone_idx];
        int j = extent.idx_to_j[zone_idx];
        int k = extent.idx_to_k[zone_idx];

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
    }); // forall3
  })); // policy
}


