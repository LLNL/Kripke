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

#ifndef KRIPKE_KERNEL_H__
#define KRIPKE_KERNEL_H__

#include <Kripke.h>
#include <Kripke/Core/DataStore.h>
#include <utility>
#include <Kripke/Arch/Scattering.h>
#include <Kripke/Arch/Source.h>
#include <Kripke/Arch/LTimes.h>
#include <Kripke/Arch/LPlusTimes.h>
#include <Kripke/Core/VarLayout.h>
#include <Kripke/Kernel.h>
#include <Kripke/Timing.h>
#include <Kripke/VarTypes.h>
#include <Kripke/Arch/LTimes.h>

namespace Kripke {

  namespace Kernel {

    using namespace Kripke;
    using namespace Kripke::Core;
    template<typename AL>
    [[clang::jit]] void LPlusTimes(Kripke::Core::DataStore &data_store){
      KRIPKE_TIMER(data_store, LPlusTimes);

      Set const &set_dir    = data_store.getVariable<Set>("Set/Direction");
      Set const &set_group  = data_store.getVariable<Set>("Set/Group");
      Set const &set_zone   = data_store.getVariable<Set>("Set/Zone");
      Set const &set_moment = data_store.getVariable<Set>("Set/Moment");

      auto &field_phi_out =   data_store.getVariable<Field_Moments>("phi_out");
      auto &field_rhs =       data_store.getVariable<Field_Flux>("rhs");
      auto &field_ell_plus =  data_store.getVariable<Field_EllPlus>("ell_plus");

      ArchLayoutV al_v = data_store.getVariable<ArchLayout>("al").al_v;

      // Loop over Subdomains
      for (Kripke::SdomId sdom_id : field_rhs.getWorkList()){
    
        using ExecPolicy = typename Kripke::Arch::Policy_LPlusTimes<AL>::ExecPolicy;

        auto sdom_al = getSdomAL(AL(), sdom_id);

    // Get dimensioning
    int num_directions = set_dir.size(sdom_id);
    int num_groups =     set_group.size(sdom_id);
    int num_moments =    set_moment.size(sdom_id);
    int num_zones =      set_zone.size(sdom_id);

    // Get views
    auto phi_out  = sdom_al.getView(field_phi_out);
    auto rhs      = sdom_al.getView(field_rhs);
    auto ell_plus = sdom_al.getView(field_ell_plus); 

    // Compute:  rhs =  ell_plus * phi_out
    RAJA::kernel<ExecPolicy>(
        camp::make_tuple(
            RAJA::TypedRangeSegment<Direction>(0, num_directions),
            RAJA::TypedRangeSegment<Moment>(0, num_moments),
            RAJA::TypedRangeSegment<Group>(0, num_groups),
            RAJA::TypedRangeSegment<Zone>(0, num_zones) ),
        KRIPKE_LAMBDA (Direction d, Moment nm, Group g, Zone z) {

            rhs(d,g,z) += ell_plus(d, nm) * phi_out(nm, g, z);

        }
    );


      }

    }


    //void LTimes(Kripke::Core::DataStore &data_store);
    
    template<typename AL>
    [[clang::jit]] void LTimes(Kripke::Core::DataStore &data_store){
        KRIPKE_TIMER(data_store, LTimes);

        using ExecPolicy = typename Kripke::Arch::Policy_LTimes<AL>::ExecPolicy;

        Set const &set_dir    = data_store.getVariable<Set>("Set/Direction");
        Set const &set_group  = data_store.getVariable<Set>("Set/Group");
        Set const &set_zone   = data_store.getVariable<Set>("Set/Zone");
        Set const &set_moment = data_store.getVariable<Set>("Set/Moment");

        auto &field_psi =       data_store.getVariable<Field_Flux>("psi");
        auto &field_phi =       data_store.getVariable<Field_Moments>("phi");
        auto &field_ell =       data_store.getVariable<Field_Ell>("ell");

        // Loop over Subdomains
        for (Kripke::SdomId sdom_id : field_psi.getWorkList()){

          auto sdom_al = getSdomAL(AL(), sdom_id);
 
          // Get dimensioning
          int num_directions = set_dir.size(sdom_id);
          int num_groups =     set_group.size(sdom_id);
          int num_moments =    set_moment.size(sdom_id);
          int num_zones =      set_zone.size(sdom_id);

          // Get pointers
          auto psi = sdom_al.getView(field_psi);
          auto phi = sdom_al.getView(field_phi);
          auto ell = sdom_al.getView(field_ell);

          // Compute:  phi =  ell * psi
          RAJA::kernel<ExecPolicy>(
              camp::make_tuple(
                  RAJA::TypedRangeSegment<Moment>(0, num_moments),
                  RAJA::TypedRangeSegment<Direction>(0, num_directions),
                  RAJA::TypedRangeSegment<Group>(0, num_groups),
                  RAJA::TypedRangeSegment<Zone>(0, num_zones) ),
              KRIPKE_LAMBDA (Moment nm, Direction d, Group g, Zone z) {

                 phi(nm,g,z) += ell(nm, d) * psi(d, g, z);

              }
          );

        }

    }

    double population(Kripke::Core::DataStore &data_store);

    template<typename AL>
    [[clang::jit]] void scattering(Kripke::Core::DataStore &data_store){
  KRIPKE_TIMER(data_store, Scattering);

  ArchLayoutV al_v = data_store.getVariable<ArchLayout>("al").al_v;

  auto &pspace = data_store.getVariable<Kripke::Core::PartitionSpace>("pspace");

  auto &set_group  = data_store.getVariable<Kripke::Core::Set>("Set/Group");
  auto &set_moment = data_store.getVariable<Kripke::Core::Set>("Set/Moment");
  auto &set_zone   = data_store.getVariable<Kripke::Core::Set>("Set/Zone");

  auto &field_phi     = data_store.getVariable<Kripke::Field_Moments>("phi");
  auto &field_phi_out = data_store.getVariable<Kripke::Field_Moments>("phi_out");
  auto &field_sigs    = data_store.getVariable<Field_SigmaS>("data/sigs");

  auto &field_zone_to_mixelem     = data_store.getVariable<Field_Zone2MixElem>("zone_to_mixelem");
  auto &field_zone_to_num_mixelem = data_store.getVariable<Field_Zone2Int>("zone_to_num_mixelem");
  auto &field_mixed_to_material = data_store.getVariable<Field_MixElem2Material>("mixelem_to_material");
  auto &field_mixed_to_fraction = data_store.getVariable<Field_MixElem2Double>("mixelem_to_fraction");

  auto &field_moment_to_legendre = data_store.getVariable<Field_Moment2Legendre>("moment_to_legendre");


  // Loop over subdomains and compute scattering source
  for(auto sdom_src : field_phi.getWorkList()){
    for(auto sdom_dst : field_phi_out.getWorkList()){

      // Only work on subdomains where src and dst are on the same R subdomain
      size_t r_src = pspace.subdomainToSpace(SPACE_R, sdom_src);
      size_t r_dst = pspace.subdomainToSpace(SPACE_R, sdom_dst);
      if(r_src != r_dst){
        continue;
      }


    using ExecPolicy = typename Kripke::Arch::Policy_Scattering<AL>::ExecPolicy;

    auto sdom_al = getSdomAL(AL(), sdom_src);

    // Get glower for src and dst ranges (to index into sigma_s)
    int glower_src = set_group.lower(sdom_src);
    int glower_dst = set_group.lower(sdom_dst);


    // get material mix information
    auto moment_to_legendre = sdom_al.getView(field_moment_to_legendre);

    auto phi     = sdom_al.getView(field_phi);
    auto phi_out = sdom_al.getView(field_phi_out, sdom_dst);
    auto sigs    = sdom_al.getView(field_sigs);
    
    auto zone_to_mixelem     = sdom_al.getView(field_zone_to_mixelem);
    auto zone_to_num_mixelem = sdom_al.getView(field_zone_to_num_mixelem);
    auto mixelem_to_material = sdom_al.getView(field_mixed_to_material);
    auto mixelem_to_fraction = sdom_al.getView(field_mixed_to_fraction);
    
    // grab dimensions
    int num_zones =      set_zone.size(sdom_src);
    int num_groups_src = set_group.size(sdom_src);
    int num_groups_dst = set_group.size(sdom_dst);
    int num_moments =    set_moment.size(sdom_dst);

    RAJA::kernel<ExecPolicy>(
        camp::make_tuple(
            RAJA::TypedRangeSegment<Moment>(0, num_moments),
            RAJA::TypedRangeSegment<Group>(0, num_groups_dst),
            RAJA::TypedRangeSegment<Group>(0, num_groups_src),
            RAJA::TypedRangeSegment<Zone>(0, num_zones) ),
        KRIPKE_LAMBDA (Moment nm, Group g, Group gp, Zone z) {

            // map nm to n
            Legendre n = moment_to_legendre(nm);

            GlobalGroup global_g{*g+glower_dst};
            GlobalGroup global_gp{*gp+glower_src};

            MixElem mix_start = zone_to_mixelem(z);
            MixElem mix_stop = mix_start + zone_to_num_mixelem(z);

            double sigs_z = 0.0;
            for(MixElem mix = mix_start;mix < mix_stop;++ mix){
              Material mat = mixelem_to_material(mix);
              double fraction = mixelem_to_fraction(mix);

              sigs_z += sigs(mat, n, global_g, global_gp) * fraction;
            }
            phi_out(nm, g, z) += sigs_z * phi(nm, gp, z);
        }
    );



    }

  }


    }

  template<typename AL>
  [[clang::jit]]  void source(Kripke::Core::DataStore &data_store){
  KRIPKE_TIMER(data_store, Source);

  ArchLayoutV al_v = data_store.getVariable<ArchLayout>("al").al_v;

  auto &set_group   = data_store.getVariable<Set>("Set/Group");
  auto &set_mixelem = data_store.getVariable<Set>("Set/MixElem");

  auto &field_phi_out = data_store.getVariable<Kripke::Field_Moments>("phi_out");

  auto &field_mixed_to_zone     = data_store.getVariable<Field_MixElem2Zone>("mixelem_to_zone");
  auto &field_mixed_to_material = data_store.getVariable<Field_MixElem2Material>("mixelem_to_material");
  auto &field_mixed_to_fraction = data_store.getVariable<Field_MixElem2Double>("mixelem_to_fraction");

  double source_strength = 1.0;


  // Loop over zoneset subdomains
  for(auto sdom_id : field_phi_out.getWorkList()){


    using ExecPolicy = typename Kripke::Arch::Policy_Source<AL>::ExecPolicy;

    auto sdom_al = getSdomAL(AL(), sdom_id);

    // Source term is isotropic
    Moment nm{0};

    auto phi_out = sdom_al.getView(field_phi_out);

    auto mixelem_to_zone     = sdom_al.getView(field_mixed_to_zone);
    auto mixelem_to_material = sdom_al.getView(field_mixed_to_material);
    auto mixelem_to_fraction = sdom_al.getView(field_mixed_to_fraction);

    int num_mixed  = set_mixelem.size(sdom_id);
    int num_groups = set_group.size(sdom_id);


    // Compute:  phi =  ell * psi
    RAJA::kernel<ExecPolicy>(
        camp::make_tuple(
            RAJA::TypedRangeSegment<Group>(0, num_groups),
            RAJA::TypedRangeSegment<MixElem>(0, num_mixed) ),
        KRIPKE_LAMBDA (Group g, MixElem mix) {

            Material material = mixelem_to_material(mix);

            if(material == 2){
              Zone z = mixelem_to_zone(mix);
              double fraction = mixelem_to_fraction(mix);

              phi_out(nm, g, z) += source_strength * fraction;
            }

        }
    );


  }


    }


    void sweepSubdomain(Kripke::Core::DataStore &data_store, Kripke::SdomId sdom_id);


    template<typename FieldType>
    RAJA_INLINE
    void kConst(FieldType &field, Kripke::SdomId sdom_id, typename FieldType::ElementType value){
      auto view1d = field.getView1d(sdom_id);
      int num_elem = field.size(sdom_id);
      RAJA::forall<RAJA::loop_exec>(
        RAJA::RangeSegment(0, num_elem),
        [=](RAJA::Index_type i){
			 	  view1d(i) = value;
      });
    }

    template<typename FieldType>
    RAJA_INLINE
    void kConst(FieldType &field, typename FieldType::ElementType value){
      for(Kripke::SdomId sdom_id : field.getWorkList()){
        kConst(field, sdom_id, value);
      }
    }




    template<typename FieldType>
    RAJA_INLINE
    void kCopy(FieldType &field_dst, Kripke::SdomId sdom_id_dst,
               FieldType &field_src, Kripke::SdomId sdom_id_src){
      auto view_src = field_src.getView1d(sdom_id_src);
      auto view_dst = field_dst.getView1d(sdom_id_dst);
      int num_elem = field_src.size(sdom_id_src);

      RAJA::forall<RAJA::loop_exec>(
        RAJA::RangeSegment(0, num_elem),
        [=](RAJA::Index_type i){
          view_src(i) = view_dst(i);
      });
    }

    template<typename FieldType>
    RAJA_INLINE
    void kCopy(FieldType &field_dst, FieldType &field_src){
      for(Kripke::SdomId sdom_id : field_dst.getWorkList()){
        kCopy(field_dst, sdom_id, field_src, sdom_id);
      }
    }

  }
}

#endif
