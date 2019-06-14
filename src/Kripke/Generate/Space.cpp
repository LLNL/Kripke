//
// Copyright (c) 2014-19, Lawrence Livermore National Security, LLC
// and Kripke project contributors. See the COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
//

#include <Kripke/Generate.h>

#include <Kripke/Core/Comm.h>
#include <Kripke/Core/Field.h>
#include <Kripke/Kernel.h>
#include <Kripke/Core/PartitionSpace.h>
#include <Kripke/Core/Set.h>
#include <Kripke/Timing.h>
#include <Kripke/VarTypes.h>

using namespace Kripke;
using namespace Kripke::Core;

namespace {

struct ZoneMixture {
  double fraction[3];

  size_t numMixed() const {

    return ( fraction[0] > 0.0 ? 1 : 0 ) +
           ( fraction[1] > 0.0 ? 1 : 0 ) +
           ( fraction[2] > 0.0 ? 1 : 0 );

  }
};

}


void Kripke::Generate::generateSpace(Kripke::Core::DataStore &data_store,
    InputVariables const &input_vars)
{
  PartitionSpace &pspace = data_store.getVariable<PartitionSpace>("pspace");

  ArchLayoutV al_v = data_store.getVariable<ArchLayout>("al").al_v;

  // Create set for X mesh
  size_t nx_per_sdom = input_vars.nx /
                       pspace.getGlobalNumSubdomains(SPACE_RX);

  KRIPKE_ASSERT(nx_per_sdom * pspace.getGlobalNumSubdomains(SPACE_RX) == (size_t)input_vars.nx,
      "Number of zones in X must evenly divide into the number of subdomains\n");

  std::vector<size_t> local_nx(pspace.getNumSubdomains(SPACE_RX),
                               nx_per_sdom);

  auto &set_zonei = data_store.newVariable<RangeSet>(
      "Set/ZoneI", pspace, SPACE_RX, local_nx);




  // Create set for Y mesh
  size_t ny_per_sdom = input_vars.ny /
                       pspace.getGlobalNumSubdomains(SPACE_RY);

  KRIPKE_ASSERT(ny_per_sdom * pspace.getGlobalNumSubdomains(SPACE_RY) == (size_t)input_vars.ny,
        "Number of zones in Y must evenly divide into the number of subdomains\n");

  std::vector<size_t> local_ny(pspace.getNumSubdomains(SPACE_RY),
                               ny_per_sdom);

  auto &set_zonej = data_store.newVariable<RangeSet>(
        "Set/ZoneJ", pspace, SPACE_RY, local_ny);



  // Create set for Z mesh
  size_t nz_per_sdom = input_vars.nz /
                       pspace.getGlobalNumSubdomains(SPACE_RZ);

  KRIPKE_ASSERT(nz_per_sdom * pspace.getGlobalNumSubdomains(SPACE_RZ) == (size_t)input_vars.nz,
          "Number of zones in Z must evenly divide into the number of subdomains\n");

  std::vector<size_t> local_nz(pspace.getNumSubdomains(SPACE_RZ),
                               nz_per_sdom);

  auto &set_zonek = data_store.newVariable<RangeSet>(
          "Set/ZoneK", pspace, SPACE_RZ, local_nz);



  // Create a total set of zones in the problem
  auto &set_zone = data_store.newVariable<ProductSet<3>>("Set/Zone",
      pspace, SPACE_R, set_zonei, set_zonej, set_zonek);


  // Create a 1d linearized set of zones
  auto &set_zone_linear = data_store.newVariable<ProductSet<1>>("Set/ZoneLinear",
      pspace, SPACE_R, set_zone);

  // Create a set of the number of materials
  data_store.newVariable<GlobalRangeSet>("Set/Material", pspace, 3);



  /* Set grid deltas for a uniform mesh (full-space, no reflecting BC's)
   * x:   -60.0  to  60.0
   * y:  -100.0  to 100.0
   * z:   -60.0  to  60.0
   */
  double const x_min = -60.0;
  double const x_max = 60.0;

  double const y_min = -100.0;
  double const y_max = 100.0;

  double const z_min = -60.0;
  double const z_max = 60.0;


  auto &field_dx = createField<Field_ZoneI2Double>(data_store, "dx", al_v, set_zonei);
  double dx = (x_max-x_min) / set_zonei.globalSize();
  Kripke::Kernel::kConst(field_dx, dx);

  auto &field_dy = createField<Field_ZoneJ2Double>(data_store, "dy", al_v, set_zonej);
  double dy = (y_max-y_min) / set_zonej.globalSize();
  Kripke::Kernel::kConst(field_dy, dy);

  auto &field_dz = createField<Field_ZoneK2Double>(data_store, "dz", al_v, set_zonek);
  double dz = (z_max-z_min) / set_zonek.globalSize();
  Kripke::Kernel::kConst(field_dz, dz);


  // Create a zone volume field (this is simple considering our uniform grid)
  double zone_volume = dx*dy*dz;
  auto &field_volume = createField<Field_Zone2Double>(data_store, "volume", al_v, set_zone_linear);
  Kripke::Kernel::kConst(field_volume, zone_volume);


  /*
   * Define a function describing the material region distribution in space
   */
  auto material_fcn = [](double x, double y, double z) -> Material {

    // Problem is defined for one octant, with reflecting boundaries
    // We "unreflect" it here by taking abs values
    x = std::abs(x);
    y = std::abs(y);
    z = std::abs(z);

    // Central 20x20x20 box is Region 1
    if(x <= 10.0 && y <= 10.0 && z <= 10.0){
      return Material{0};
    }

    // Leg 1 of Region 2
    if(x <= 10.0 && y <= 60.0 && z <= 10.0){
      return Material{1};
    }

    // Leg 2 of Region 2
    if(x <= 40.0 && y >= 50.0 && y <= 60.0 && z <= 10.0){
      return Material{1};
    }

    // Leg 3 of Region 2
    if(x >= 30.0 && x <= 40.0 && y >= 50.0 && y <= 60.0 && z <= 40.0){
      return Material{1};
    }

    // Leg 4 of Region 2
    if(x >= 30.0 && x <= 40.0 && y >= 50.0 && z >= 30.0 && z <= 40.0){
      return Material{1};
    }

    // Rest is filled with region 3
    return Material{2};
  };





  /*
   *  For each subdomain in space, build up dynamic arrays to hold material
   *  mixture information
   */

  // number of subsamples per spatial dimension
  int num_subsamples = input_vars.num_material_subsamples;
  double sample_vol_frac = 1.0 / (double)(num_subsamples*num_subsamples*num_subsamples);

  auto sdom_list = set_zone.getWorkList();
  std::vector<std::vector<ZoneMixture>> mix;

  for(SdomId sdom_id : sdom_list){

    double x0 = x_min + dx*set_zonei.lower(sdom_id);
    double y0 = y_min + dy*set_zonej.lower(sdom_id);
    double z0 = z_min + dz*set_zonek.lower(sdom_id);

    std::vector<ZoneMixture> sdom_mix(set_zone.size(sdom_id));
    auto zone_layout = set_zone.getLayout(sdom_id);

    // iterate over the zones, assume uniform mesh for our coordinate
    // calculations
    for (int i = 0; i < (int)set_zonei.size(sdom_id); i ++) {
      for (int j = 0; j < (int)set_zonej.size(sdom_id); j ++) {
        for (int k = 0; k < (int)set_zonek.size(sdom_id); k ++) {

          int zone = zone_layout(i,j,k);

          double xi = x0 + dx*i;
          double yi = y0 + dy*j;
          double zi = z0 + dz*k;

          // subsample probe the geometry to get our materials
          sdom_mix[zone] = {{0.0, 0.0, 0.0}}; // fraction of both materials

          for(int si = 0;si < num_subsamples;++ si){
            for(int sj = 0;sj < num_subsamples;++ sj){
              for(int sk = 0;sk < num_subsamples;++ sk){

                double x = xi + dx*(si+1)/(num_subsamples+1);
                double y = yi + dy*(sj+1)/(num_subsamples+1);
                double z = zi + dz*(sk+1)/(num_subsamples+1);

                Material mat = material_fcn(x, y, z);
                sdom_mix[zone].fraction[*mat] += sample_vol_frac;

              }
            }
          }
        }
      }
    }

    mix.push_back(sdom_mix);
  } // sdom_id


  // Go through and count number of mixed elements
  std::vector<size_t> sdom_to_num_mixed;
  for(auto &sdom_mix : mix){
    size_t n = 0;
    for(auto &z : sdom_mix){
      n += z.numMixed();
    }
    sdom_to_num_mixed.push_back(n);
  }

  // Create a new set that describes the number of mixed zones per sdom
  auto &set_mixelem = data_store.newVariable<RangeSet>(
            "Set/MixElem", pspace, SPACE_R, sdom_to_num_mixed);

  // Create fields to store mixture information
  auto &field_mixed_to_zone = createField<Field_MixElem2Zone>(
      data_store, "mixelem_to_zone", al_v, set_mixelem);

  auto &field_mixed_to_material = createField<Field_MixElem2Material>(
      data_store, "mixelem_to_material", al_v, set_mixelem);

  auto &field_mixed_to_fraction = createField<Field_MixElem2Double>(
      data_store, "mixelem_to_fraction", al_v, set_mixelem);

  auto &field_zone_to_num_mixelem = createField<Field_Zone2Int>(
      data_store, "zone_to_num_mixelem", al_v, set_zone_linear);

  auto &field_zone_to_mixelem = createField<Field_Zone2MixElem>(
      data_store, "zone_to_mixelem", al_v, set_zone_linear);


  // Populate mixture fields with our dynamic data
  RAJA::ReduceSum<RAJA::seq_reduce, double> total_volume_red[3];
  for(size_t i = 0;i < sdom_list.size();++ i){
    SdomId sdom_id = sdom_list[i];

    int num_zones = set_zone.size(sdom_id);
    int num_mixelems = set_mixelem.size(sdom_id);

    auto const &sdom_mix = mix[i];

    auto mixed_to_zone       = field_mixed_to_zone.getView(sdom_id);
    auto mixed_to_material   = field_mixed_to_material.getView(sdom_id);
    auto mixed_to_fraction   = field_mixed_to_fraction.getView(sdom_id);
    auto zone_to_num_mixelem = field_zone_to_num_mixelem.getView(sdom_id);
    auto zone_to_mixelem     = field_zone_to_mixelem.getView(sdom_id);

    RAJA::ReduceSum<RAJA::seq_reduce, MixElem> mixelem(MixElem{0});
    RAJA::forall<RAJA::seq_exec>(
      RAJA::TypedRangeSegment<Zone>(0, num_zones),
      [=](Zone z){
        ZoneMixture const &zone_mix = sdom_mix[*z];
        int num_zone_mix = 0;

        zone_to_mixelem(z) = mixelem;

        double zone_frac = 0.0;
        for(Material m{0};m < 3;++ m){
          if(zone_mix.fraction[*m] > 0.0){
            MixElem me = mixelem;

            mixed_to_zone(me) = z;
            mixed_to_material(me) = m;
            mixed_to_fraction(me) = zone_mix.fraction[*m];
            zone_frac += zone_mix.fraction[*m];
            total_volume_red[*m] += zone_mix.fraction[*m] * zone_volume;
            num_zone_mix ++;
            mixelem += MixElem{1};
          }
        }
        KRIPKE_ASSERT(zone_frac == 1.0, "Zone fraction wrong: %e", zone_frac);
        zone_to_num_mixelem(z) = num_zone_mix;
    });

    KRIPKE_ASSERT((*((MixElem)mixelem)) == num_mixelems, "Mismatch in mixture info");
  }

  // Display the total volume
  Kripke::Core::Comm default_comm;
  auto const &r_comm = pspace.getComm(SPACE_R);
  double total_volume[3];
  total_volume[0] = total_volume_red[0];
  total_volume[1] = total_volume_red[1];
  total_volume[2] = total_volume_red[2];
  r_comm.allReduceSumDouble(total_volume, 3);
  if(default_comm.rank() == 0){
    printf("\n  Material Volumes=[%e, %e, %e]\n", total_volume[0],
      total_volume[1],total_volume[2]);
  }


  // Allocate storage for our zonal total cross-section
  auto &set_group = data_store.getVariable<Set>("Set/Group");
  auto &set_sigt_zonal = data_store.newVariable<ProductSet<2>>(
      "Set/SigmaTZonal", pspace, SPACE_PR, set_group, set_zone);
  auto &field_sigt = createField<Field_SigmaTZonal>(
      data_store, "sigt_zonal", al_v, set_sigt_zonal);
  Kripke::Kernel::kConst(field_sigt, 0.0);

  for(SdomId sdom_id : field_sigt.getWorkList()){

    auto mixelem_to_zone     = field_mixed_to_zone.getView(sdom_id);
    auto mixelem_to_material = field_mixed_to_material.getView(sdom_id);
    auto mixelem_to_fraction = field_mixed_to_fraction.getView(sdom_id);
    auto sigt = field_sigt.getView(sdom_id);

    int num_groups  = set_group.size(sdom_id);
    int num_mixelem = set_mixelem.size(sdom_id);

    for(Group g{0};g < num_groups;++ g){

      RAJA::forall<RAJA::seq_exec>(
        RAJA::TypedRangeSegment<MixElem>(0, num_mixelem),
        [=](MixElem mixelem){
          Zone z       = mixelem_to_zone(mixelem);
          Material mat = mixelem_to_material(mixelem);

          sigt(g, z) += mixelem_to_fraction(mixelem) * input_vars.sigt[*mat];
      });
    }

  }

}




