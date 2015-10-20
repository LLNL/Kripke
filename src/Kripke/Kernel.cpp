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
void Kernel::LTimes(Grid_Data *domain) { 
  
  policyScope(nesting_order, [&](auto nest_tag){
    typedef decltype(nest_tag) nest_type;
    typedef LayoutPolicy<nest_type> LAYOUT;
    typedef ViewPolicy<LAYOUT> VIEW;

    // Zero Phi
    forallZoneSets(domain, [&](int zs, int sdom_id, Subdomain &sdom){
      sdom.phi->clear(0.0);
    });

    // Loop over Subdomains
    forallSubdomains(domain, [&](int sdom_id, Subdomain &sdom){

      // Get dimensioning
      int group0 = sdom.group0;

      // Get pointers
      typename VIEW::TPsi psi(sdom.psi->ptr(), domain, sdom_id);
      typename VIEW::TPhi phi(sdom.phi->ptr(), domain, sdom_id);
      typename VIEW::TEll ell(sdom.ell->ptr(), domain, sdom_id);

      forall4T<LTimesPolicy<nest_type>, IMoment, IDirection, IGroup, IZone>(
        domain, sdom_id, 
        [&](IMoment nm, IDirection d, IGroup g, IZone z){
  
          IGlobalGroup g_global( (*g) + group0);
          
          phi(nm, g_global, z) += ell(d, nm) * psi(d, g, z);
      });

    }); // sdom

  }); // policy
}

#include<Kripke/Kernel/LPlusTimesPolicy.h>
void Kernel::LPlusTimes(Grid_Data *domain) {
  policyScope(nesting_order, [&](auto nest_tag){
    typedef decltype(nest_tag) nest_type;
    typedef LayoutPolicy<nest_type> LAYOUT;
    typedef ViewPolicy<LAYOUT> VIEW;
    // Loop over Subdomains
    forallSubdomains(domain, [&](int sdom_id, Subdomain &sdom){
      sdom.rhs->clear(0.0);
    });

    forallSubdomains(domain, [&](int sdom_id, Subdomain &sdom){

      // Get dimensioning
      int group0 = sdom.group0;

      // Get pointers
      typename VIEW::TPsi     rhs(sdom.rhs->ptr(), domain, sdom_id);
      typename VIEW::TPhi     phi_out(sdom.phi_out->ptr(), domain, sdom_id);
      typename VIEW::TEllPlus ell_plus(sdom.ell_plus->ptr(), domain, sdom_id);
      
      forall4T<LPlusTimesPolicy<nest_type>, IMoment, IDirection, IGroup, IZone>(
        domain, sdom_id, 
        [&](IMoment nm, IDirection d, IGroup g, IZone z){
  
          IGlobalGroup g_global( (*g) + group0);
  
          rhs(d, g, z) += ell_plus(d, nm) * phi_out(nm, g_global, z);  
      });

    }); // sdom

  }); // policy
}


/**
  Compute scattering source term phi_out from flux moments in phi.
  phi_out(gp,z,nm) = sum_g { sigs(g, n, gp) * phi(g,z,nm) }
*/
#include<Kripke/Kernel/ScatteringPolicy.h>
void Kernel::scattering(Grid_Data *domain){
  
  policyScope(nesting_order, [&](auto nest_tag){
    typedef decltype(nest_tag) nest_type;
    typedef LayoutPolicy<nest_type> LAYOUT;
    typedef ViewPolicy<LAYOUT> VIEW;

    // Zero out source terms
    forallZoneSets(domain, [&](int zs, int sdom_id, Subdomain &sdom){
      sdom.phi_out->clear(0.0);
    });

    // Loop over zoneset subdomains
    forallZoneSets(domain, [&](int zs, int sdom_id, Subdomain &sdom){

      typename VIEW::TPhi     phi(sdom.phi->ptr(), domain, sdom_id);
      typename VIEW::TPhi     phi_out(sdom.phi_out->ptr(), domain, sdom_id);
      typename VIEW::TSigS    sigs(domain->sigs->ptr(), domain, sdom_id);

      View1d<const int, PERM_I>    const mixed_to_zones(&sdom.mixed_to_zones[0], 1);
      View1d<const int, PERM_I>    const mixed_material(&sdom.mixed_material[0], 1);
      View1d<const double, PERM_I> const mixed_fraction(&sdom.mixed_fraction[0], 1);
      View1d<const int, PERM_I>    const moment_to_coeff(&domain->moment_to_coeff[0], 1);
      
      forall4T<ScatteringPolicy<nest_type>, IMoment, IGlobalGroup, IGlobalGroup, IMix>(
        domain, sdom_id,
        [&](IMoment nm, IGlobalGroup g, IGlobalGroup gp, IMix mix){
        
          ILegendre n(moment_to_coeff(*nm));
          IZone zone(mixed_to_zones(*mix));
          IMaterial material(mixed_material(*mix));
          double fraction = mixed_fraction(*mix);

          phi_out(nm, gp, zone) += 
            sigs(n, g, gp, material) * phi(nm, g, zone) * fraction;                     
                             
        });  // forall
          
    }); // zonesets

  }); // policy
}

  
/**
 * Add an isotropic source, with flux of 1, to every zone with Region 1
 * (or material 0).
 *
 * Since it's isotropic, we're just adding this to nm=0.
 */
#include<Kripke/Kernel/SourcePolicy.h>
void Kernel::source(Grid_Data *grid_data){

  policyScope(nesting_order, [&](auto nest_tag){
    typedef decltype(nest_tag) nest_type;
    typedef LayoutPolicy<nest_type> LAYOUT;
    typedef ViewPolicy<LAYOUT> VIEW;

    // Loop over zoneset subdomains
    forallZoneSets(grid_data, [&](int zs, int sdom_id, Subdomain &sdom){
      typename VIEW::TPhi     phi_out(sdom.phi_out->ptr(), grid_data, sdom_id);
      View1d<const int,    PERM_I> const mixed_to_zones(&sdom.mixed_to_zones[0], 1);
      View1d<const int,    PERM_I> const mixed_material(&sdom.mixed_material[0], 1);
      View1d<const double, PERM_I> const mixed_fraction(&sdom.mixed_fraction[0], 1);

      forall2T<SourcePolicy<nest_type>, IGlobalGroup, IMix >(
        grid_data, sdom_id,
        [&](IGlobalGroup g, IMix mix){
          int zone = mixed_to_zones(*mix);
          int material = mixed_material(*mix);
          double fraction = mixed_fraction(*mix);

          if(material == 0){
            phi_out(IMoment(0), g, IZone(zone)) += 1.0 * fraction;
          }
      }); // forall

    }); // zonesets

  }); // policy
}


#include<Kripke/Kernel/SweepPolicy.h>
void Kernel::sweep(Grid_Data *domain, int sdom_id) {
  policyScope(nesting_order, [&](auto nest_tag){
    typedef decltype(nest_tag) nest_type;
    typedef LayoutPolicy<nest_type> LAYOUT;
    typedef ViewPolicy<LAYOUT> VIEW;

    Subdomain *sdom = &domain->subdomains[sdom_id];

    int num_directions = sdom->num_directions;
    int num_groups = sdom->num_groups;
    int num_zones = sdom->num_zones;

    int local_imax = sdom->nzones[0];
    int local_jmax = sdom->nzones[1];
    int local_kmax = sdom->nzones[2];

    int num_z_i = local_jmax * local_kmax;
    int num_z_j = local_imax * local_kmax;
    int num_z_k = local_imax * local_jmax;

    typename VIEW::TDirections direction(sdom->directions, domain, sdom_id);

    typename VIEW::TPsi     rhs(sdom->rhs->ptr(), domain, sdom_id);
    typename VIEW::TPsi     psi(sdom->psi->ptr(), domain, sdom_id);
    typename VIEW::TSigT    sigt(sdom->sigt->ptr(), domain, sdom_id);

    typename VIEW::Tdx      dx(&sdom->deltas[0][0], domain, sdom_id);
    typename VIEW::Tdy      dy(&sdom->deltas[1][0], domain, sdom_id);
    typename VIEW::Tdz      dz(&sdom->deltas[2][0], domain, sdom_id);
    
    typename LAYOUT::TZone zone_layout(domain, sdom_id);
    
    typename VIEW::TFaceI face_lf(sdom->plane_data[0]->ptr(), domain, sdom_id);
    typename VIEW::TFaceJ face_fr(sdom->plane_data[1]->ptr(), domain, sdom_id);
    typename VIEW::TFaceK face_bo(sdom->plane_data[2]->ptr(), domain, sdom_id);

    // All directions have same id,jd,kd, since these are all one Direction Set
    // So pull that information out now
    Grid_Sweep_Block const &extent = sdom->sweep_block;    
    typename VIEW::TIdxToI  idx_to_i((IZoneI*)&extent.idx_to_i[0], domain, sdom_id);
    typename VIEW::TIdxToJ  idx_to_j((IZoneJ*)&extent.idx_to_j[0], domain, sdom_id);
    typename VIEW::TIdxToK  idx_to_k((IZoneK*)&extent.idx_to_k[0], domain, sdom_id);

    forall3T<SweepPolicy<nest_type>, IDirection, IGroup, IZoneIdx>(
      IDirection::range(domain, sdom_id),
      IGroup::range(domain, sdom_id),
      extent.indexset_sweep,
      [&](IDirection d, IGroup g, IZoneIdx zone_idx){
        
        IZoneI i = idx_to_i(zone_idx);
        IZoneJ j = idx_to_j(zone_idx);
        IZoneK k = idx_to_k(zone_idx);
        
        double const xcos_dxi = 2.0 * direction(d).xcos / dx(i+1); 
        double const ycos_dyj = 2.0 * direction(d).ycos / dy(j+1); 
        double const zcos_dzk = 2.0 * direction(d).zcos / dz(k+1); 

        IZone z = zone_layout(i,j,k);
        
        /* Calculate new zonal flux */
        double const psi_d_g_z = (
              rhs(d,g,z)
            + face_lf(d,g,j,k) * xcos_dxi
            + face_fr(d,g,i,k) * ycos_dyj
            + face_bo(d,g,i,j) * zcos_dzk)
            / (xcos_dxi + ycos_dyj + zcos_dzk + sigt(g,z) );

        psi(d,g,z) = psi_d_g_z;

        /* Apply diamond-difference relationships */
        face_lf(d,g,j,k) = 2.0 * psi_d_g_z - face_lf(d,g,j,k);
        face_fr(d,g,i,k) = 2.0 * psi_d_g_z - face_fr(d,g,i,k);
        face_bo(d,g,i,j) = 2.0 * psi_d_g_z - face_bo(d,g,i,j);
    }); // forall3
  }); // policy  
}

