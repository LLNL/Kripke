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

#include <Kripke/Generate.h>

#include <Kripke/Comm.h>
#include <Kripke/Field.h>
#include <Kripke/Kernel.h>
#include <Kripke/Quadrature.h>
#include <Kripke/PartitionSpace.h>
#include <Kripke/Set.h>
#include <Kripke/Timing.h>
#include <Kripke/VarTypes.h>

using namespace Kripke;


static void initializeDecomp(Kripke::DataStore &data_store,
    InputVariables const &input_vars)
{
  // Create a "Comm World"
  auto &comm = data_store.newVariable<Comm>("comm");

  // Create our partitioning over MPI
  auto &pspace = data_store.newVariable<PartitionSpace>("pspace",
      comm,
      1,
      1,
      input_vars.npx,
      input_vars.npy,
      input_vars.npz);

  // Create our local partition over subdomains
  pspace.setup_createSubdomains(
      input_vars.num_groupsets,
      input_vars.num_dirsets,
      input_vars.num_zonesets_dim[0],
      input_vars.num_zonesets_dim[1],
      input_vars.num_zonesets_dim[2]);

  // Create utility Sets and Fields that describe our global subdomain layout
  pspace.createSubdomainData(data_store);

  pspace.print();

}





static void initializeEnergy(Kripke::DataStore &data_store,
    InputVariables const &input_vars)
{

  PartitionSpace &pspace = data_store.getVariable<PartitionSpace>("pspace");

  // Create sets for energy discretization
  size_t ngrp_per_sdom = input_vars.num_groups /
                         pspace.getGlobalNumSubdomains(SPACE_P);

  std::vector<size_t> local_grps(pspace.getNumSubdomains(SPACE_P),
                                 ngrp_per_sdom);

  RangeSet *grp_set = new RangeSet(pspace, SPACE_P, local_grps);
  data_store.addVariable("Set/Group", grp_set);

  GlobalRangeSet *global_grp_set = new GlobalRangeSet(pspace, *grp_set);
  data_store.addVariable("Set/GlobalGroup", global_grp_set);


}



namespace {
  double FactFcn(int n)
  {
    double fact = 1.0;
    for(int i = n;i > 0 ;--i){
      fact *= (double)i;
    }
    return(fact);
  }

  inline double PnmFcn(int n, int m, double x)
  {
    /*-----------------------------------------------------------------
     * It is assumed that 0 <= m <= n and that abs(x) <= 1.0.
     * No error checking is done, however.
     *---------------------------------------------------------------*/
    double fact, pnn=0, pmm, pmmp1, somx2;

    int i, nn;

    if(std::abs(x) > 1.0){
      KRIPKE_ABORT("Bad input to ardra_PnmFcn: abs(x) > 1.0, x = %e\n", x);
    }
    else if((x > 1.0) && (x <= 1.0)){
      x = 1.0;
    }
    else if((-1.0 <= x ) && (x < -1.0)){
      x = -1.0;
    }

    pmm=1.0;
    if(m > 0){
      somx2=sqrt((1.0-x)*(1.0+x));
      fact=1.0;
      for(i=1; i<=m; i++){
        pmm *= -fact*somx2;
        fact += 2.0;
      }
    }
    if(n == m){
      return(pmm);
    }
    else {
      pmmp1=x*(2*m+1)*pmm;
      if(n == (m+1)){
        return(pmmp1);
      }
      else {
        for(nn=m+2; nn<=n; nn++){
          pnn=(x*(2*nn-1)*pmmp1-(nn+m-1)*pmm)/(nn-m);
          pmm=pmmp1;
          pmmp1=pnn;
        }
        return(pnn);
      }
    }
  }

  inline double YnmFcn(int n, int m, double mu, double eta, double xi)
  {
    double fac1, fac2, anm, ynm, pnm, dm0, taum, tmp, phi, phi_tmp;
    double floor=1.e-20;
    int nn, mm;

    /* Calculate the correct phi for omega=(mu,eta,xi) */
    tmp = fabs(eta/(mu+floor));
    phi_tmp = atan(tmp);
    if( (mu>0) && (eta>0) ){
      phi = phi_tmp;
    }
    else if( (mu<0) && (eta>0) ){
      phi = M_PI - fabs(phi_tmp);
    }
    else if( (mu<0) && (eta<0) ){
      phi = M_PI + fabs(phi_tmp);
    }
    else {
      phi = 2.0*M_PI - fabs(phi_tmp);
    }

    /* Begin evaluation of Ynm(omega) */
    nn = n - std::abs(m);
    fac1 = (double) FactFcn(nn);
    nn = n + std::abs(m);
    fac2 = (double) FactFcn(nn);
    mm = std::abs(m);
    pnm = PnmFcn(n, mm, xi);
    tmp = ((double) m)*phi;
    if(m >= 0){
      taum = cos(tmp);
    }
    else {taum = sin(-tmp); }
    if(m == 0){
      dm0 = 1.0;
    }
    else {dm0 = 0.0; }

    tmp = ((2*n+1)*fac1)/(2.0*(1.0+dm0)*M_PI*fac2);
    anm = sqrt( tmp );
    ynm = anm*pnm*taum;
    return(ynm);
  }
}






static void initializeDirections(Kripke::DataStore &data_store,
    InputVariables const &input_vars)
{

  PartitionSpace &pspace = data_store.getVariable<PartitionSpace>("pspace");

  // Create sets for angular discretization
  size_t ndir_per_sdom = input_vars.num_directions /
                        pspace.getGlobalNumSubdomains(SPACE_Q);

  std::vector<size_t> local_dirs(pspace.getNumSubdomains(SPACE_Q),
                                 ndir_per_sdom);

  RangeSet *dir_set = new RangeSet(pspace, SPACE_Q, local_dirs);

  data_store.addVariable("Set/Direction", dir_set);


  size_t legendre_order = input_vars.legendre_order;
  size_t num_moments = (legendre_order+1)*(legendre_order+1);

  GlobalRangeSet *moment_set = new GlobalRangeSet(pspace, num_moments);
  data_store.addVariable("Set/Moment", moment_set);

  GlobalRangeSet *legendre_set = new GlobalRangeSet(pspace, legendre_order+1);
  data_store.addVariable("Set/Legendre", legendre_set);

  // Create a mapping from moments to their Legendre scattering coefficients
  auto &field_moment_to_legendre = data_store.newVariable<Field_Moment2Legendre>(
      "moment_to_legendre", *moment_set);

  // create a global mapping... just easier this way
  std::vector<Legendre> moment_list(moment_set->globalSize());
  int nm = 0;
  for(int n = 0;n < (int)legendre_order+1;++ n){
    for(int m = -n;m <= n; ++ m){
      moment_list[nm] = Legendre{n};
      ++ nm;
    }
  }
  KRIPKE_ASSERT(nm == (int)moment_set->globalSize());

  // fill in the global
  printf("num sdom=%d\n", (int)field_moment_to_legendre.getWorkList().size());
  for(SdomId sdom_id : field_moment_to_legendre.getWorkList()){
    auto moment_to_legendre = field_moment_to_legendre.getView(sdom_id);

    for(Moment nm{0};nm < moment_set->size(sdom_id);++ nm){
      moment_to_legendre(nm) = moment_list[(*nm) + moment_set->lower(sdom_id)];
    }
  }


  // Create the quadrature set
  auto quadrature_points = Kripke::createQuadratureSet(input_vars);

  // Create and populate fields for quadrature set data
  auto &field_xcos = data_store.newVariable<Field_Direction2Double>("quadrature/xcos", *dir_set);
  auto &field_ycos = data_store.newVariable<Field_Direction2Double>("quadrature/ycos", *dir_set);
  auto &field_zcos = data_store.newVariable<Field_Direction2Double>("quadrature/zcos", *dir_set);
  auto &field_w = data_store.newVariable<Field_Direction2Double>("quadrature/w", *dir_set);
  auto &field_id = data_store.newVariable<Field_Direction2Int>("quadrature/id", *dir_set);
  auto &field_jd = data_store.newVariable<Field_Direction2Int>("quadrature/jd", *dir_set);
  auto &field_kd = data_store.newVariable<Field_Direction2Int>("quadrature/kd", *dir_set);
  auto &field_octant = data_store.newVariable<Field_Direction2Int>("quadrature/octant", *dir_set);

  for(SdomId sdom_id : field_xcos.getWorkList()){
    int num_directions = dir_set->size(sdom_id);
    int direction_lower = dir_set->lower(sdom_id);

    auto xcos = field_xcos.getView(sdom_id);
    auto ycos = field_ycos.getView(sdom_id);
    auto zcos = field_zcos.getView(sdom_id);
    auto w = field_w.getView(sdom_id);
    auto id = field_id.getView(sdom_id);
    auto jd = field_jd.getView(sdom_id);
    auto kd = field_kd.getView(sdom_id);
    auto octant = field_octant.getView(sdom_id);

    for(Direction d{0};d < num_directions;++ d){
      QuadraturePoint const &point_d = quadrature_points[(*d)+direction_lower];
      xcos(d) = point_d.xcos;
      ycos(d) = point_d.ycos;
      zcos(d) = point_d.zcos;
      w(d) = point_d.w;
      id(d) = point_d.id;
      jd(d) = point_d.jd;
      kd(d) = point_d.kd;
      octant(d) = point_d.octant;
    }
  }


  // Create a set to describe the L and L+ matrices
  auto &set_ell = data_store.newVariable<ProductSet<2>>("Set/Ell",
      pspace, SPACE_Q, *moment_set, *dir_set);

  // Allocate and initialize the L and L+ matrices
  auto &field_ell = data_store.newVariable<Field_Ell>("ell", set_ell);
  auto &field_ell_plus = data_store.newVariable<Field_Ell>("ell_plus", set_ell);

  for(SdomId sdom_id : field_xcos.getWorkList()){
    auto ell = field_ell.getView(sdom_id);
    auto ell_plus = field_ell_plus.getView(sdom_id);

    int num_directions = dir_set->size(sdom_id);
    int direction_lower = dir_set->lower(sdom_id);

    double SQRT4PI = std::sqrt(4*M_PI);
    Moment nm{0};
    for(int n=0; n < (int)legendre_order+1; n++){
      for(int m=-n; m<=n; m++){
        for(Direction d{0}; d<num_directions; d++){

          QuadraturePoint const &point_d = quadrature_points[(*d)+direction_lower];
          // Get quadrature point info
          double xcos = (point_d.id)*(point_d.xcos);
          double ycos = (point_d.jd)*(point_d.ycos);
          double zcos = (point_d.kd)*(point_d.zcos);
          double w =  point_d.w;

          double ynm = YnmFcn(n, m, xcos, ycos, zcos);

          // Compute element of L and L+
          ell(nm, d) = w*ynm/SQRT4PI;
          ell_plus(nm,d) = ynm*SQRT4PI;
        }
        nm ++;
      }
    }
  }
}



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


static void initializeSpace(Kripke::DataStore &data_store,
    InputVariables const &input_vars)
{
  PartitionSpace &pspace = data_store.getVariable<PartitionSpace>("pspace");



  // Create set for X mesh
  size_t nx_per_sdom = input_vars.nx /
                       pspace.getGlobalNumSubdomains(SPACE_RX);

  std::vector<size_t> local_nx(pspace.getNumSubdomains(SPACE_RX),
                               nx_per_sdom);

  auto &set_zonei = data_store.newVariable<RangeSet>(
      "Set/ZoneI", pspace, SPACE_RX, local_nx);




  // Create set for Y mesh
  size_t ny_per_sdom = input_vars.ny /
                       pspace.getGlobalNumSubdomains(SPACE_RY);

  std::vector<size_t> local_ny(pspace.getNumSubdomains(SPACE_RY),
                               ny_per_sdom);

  auto &set_zonej = data_store.newVariable<RangeSet>(
        "Set/ZoneJ", pspace, SPACE_RY, local_ny);



  // Create set for Z mesh
  size_t nz_per_sdom = input_vars.nz /
                       pspace.getGlobalNumSubdomains(SPACE_RZ);

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
/*
  double const x_min = 0;
  double const x_max = 1;

  double const y_min = 0;
  double const y_max = 1;

  double const z_min = 0;
  double const z_max = 1; */
  auto &field_dx = data_store.newVariable<Field_ZoneI2Double>("dx", set_zonei);
  double dx = (x_max-x_min) / set_zonei.globalSize();
  Kripke::Kernel::kConst(field_dx, dx);

  auto &field_dy = data_store.newVariable<Field_ZoneJ2Double>("dy", set_zonej);
  double dy = (y_max-y_min) / set_zonej.globalSize();
  Kripke::Kernel::kConst(field_dy, dy);

  auto &field_dz = data_store.newVariable<Field_ZoneK2Double>("dz", set_zonek);
  double dz = (z_max-z_min) / set_zonek.globalSize();
  Kripke::Kernel::kConst(field_dz, dz);


  // Create a zone volume field (this is simple considering our uniform grid)
  double zone_volume = dx*dy*dz;
  auto &field_volume = data_store.newVariable<Field_Zone2Double>("volume", set_zone_linear);
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
          sdom_mix[zone] = {0.0, 0.0, 0.0}; // fraction of both materials

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

  printf("Global size of Set/MixElem: %d\n", (int)set_mixelem.globalSize());

  // Create fields to store mixture information
  auto &field_mixed_to_zone = data_store.newVariable<Field_MixElem2Zone>(
      "mixelem_to_zone", set_mixelem);

  auto &field_mixed_to_material = data_store.newVariable<Field_MixElem2Material>(
      "mixelem_to_material", set_mixelem);

  auto &field_mixed_to_fraction = data_store.newVariable<Field_MixElem2Double>(
      "mixelem_to_fraction", set_mixelem);

  auto &field_zone_to_num_mixelem = data_store.newVariable<Field_Zone2Int>(
      "zone_to_num_mixelem", set_zone_linear);

  auto &field_zone_to_mixelem = data_store.newVariable<Field_Zone2MixElem>(
      "zone_to_mixelem", set_zone_linear);


  // Populate mixture fields with our dynamic data
  double total_volume[3] = {0.0, 0.0, 0.0};
  for(size_t i = 0;i < sdom_list.size();++ i){
    SdomId sdom_id = sdom_list[i];

    int num_zones = set_zone.size(sdom_id);
    int num_mixelems = set_mixelem.size(sdom_id);

    auto &sdom_mix = mix[i];

    auto mixed_to_zone       = field_mixed_to_zone.getView(sdom_id);
    auto mixed_to_material   = field_mixed_to_material.getView(sdom_id);
    auto mixed_to_fraction   = field_mixed_to_fraction.getView(sdom_id);
    auto zone_to_num_mixelem = field_zone_to_num_mixelem.getView(sdom_id);
    auto zone_to_mixelem     = field_zone_to_mixelem.getView(sdom_id);

    MixElem mixelem{0};
    for(Zone z{0};z < num_zones;++ z){
      ZoneMixture &zone_mix = sdom_mix[*z];
      int num_zone_mix = 0;

      zone_to_mixelem(z) = mixelem;

      double zone_frac = 0.0;
      for(Material m{0};m < 3;++ m){
        if(zone_mix.fraction[*m] > 0.0){
          mixed_to_zone(mixelem) = z;
          mixed_to_material(mixelem) = m;
          mixed_to_fraction(mixelem) = zone_mix.fraction[*m];
          zone_frac += zone_mix.fraction[*m];
          total_volume[*m] += zone_mix.fraction[*m] * zone_volume;
          num_zone_mix ++;
          mixelem ++;
        }
      }
      KRIPKE_ASSERT(zone_frac == 1.0, "Zone fraction wrong: %e", zone_frac);
      zone_to_num_mixelem(z) = num_zone_mix;
    }

    KRIPKE_ASSERT((*mixelem) == num_mixelems, "Mismatch in mixture info");
  }

  printf("Volumes: %.12e, %.12e, %.12e\n", total_volume[0],
      total_volume[1],total_volume[2]);


  // Allocate storage for our zonal total cross-section
  auto &set_group = data_store.getVariable<Set>("Set/Group");
  auto &set_sigt_zonal = data_store.newVariable<ProductSet<2>>(
      "Set/SigmaTZonal", pspace, SPACE_PR, set_group, set_zone);
  auto &field_sigt = data_store.newVariable<Field_SigmaTZonal>(
      "sigt_zonal", set_sigt_zonal);
  Kripke::Kernel::kConst(field_sigt, 0.0);

  for(SdomId sdom_id : field_sigt.getWorkList()){

    auto mixelem_to_zone     = field_mixed_to_zone.getView(sdom_id);
    auto mixelem_to_material = field_mixed_to_material.getView(sdom_id);
    auto mixelem_to_fraction = field_mixed_to_fraction.getView(sdom_id);
    auto sigt = field_sigt.getView(sdom_id);

    int num_groups  = set_group.size(sdom_id);
    int num_mixelem = set_mixelem.size(sdom_id);

    for(Group g{0};g < num_groups;++ g){
      for(MixElem mixelem{0};mixelem < num_mixelem;++ mixelem){
        Zone z       = mixelem_to_zone(mixelem);
        Material mat = mixelem_to_material(mixelem);

        sigt(g, z) += mixelem_to_fraction(mixelem) * input_vars.sigt[*mat];

      }
    }

  }

}






static void initializeData(Kripke::DataStore &data_store,
    InputVariables const &input_vars)
{

  PartitionSpace &pspace = data_store.getVariable<PartitionSpace>("pspace");


  // Create a set to span angular the flux
  Set const &dir_set   = data_store.getVariable<Set>("Set/Direction");
  Set const &group_set = data_store.getVariable<Set>("Set/Group");
  Set const &zone_set  = data_store.getVariable<Set>("Set/Zone");
  ProductSet<3> *flux_set = new ProductSet<3>(pspace, SPACE_PQR,
      dir_set, group_set, zone_set);

  data_store.addVariable("Set/Flux", flux_set);

  // Create Solution and RHS fields
  data_store.addVariable("psi", new Field_Flux(*flux_set));
  data_store.addVariable("rhs", new Field_Flux(*flux_set));




  // Create a set to span moments of the angular flux
  Set const &moment_set   = data_store.getVariable<Set>("Set/Moment");
  ProductSet<3> *fluxmoment_set = new ProductSet<3>(pspace, SPACE_PR,
        moment_set, group_set, zone_set);

  data_store.addVariable("Set/FluxMoment", fluxmoment_set);


  // Create flux moment and source moment fields
  data_store.addVariable("phi",     new Field_Moments(*fluxmoment_set));
  data_store.addVariable("phi_out", new Field_Moments(*fluxmoment_set));


  // Create "plane data" to hold face-centered values while sweeping
  Set const &zonei_set = data_store.getVariable<Set>("Set/ZoneI");
  Set const &zonej_set = data_store.getVariable<Set>("Set/ZoneJ");
  Set const &zonek_set = data_store.getVariable<Set>("Set/ZoneK");
  Set const &iplane_set = data_store.newVariable<ProductSet<4>>("Set/IPlane", pspace, SPACE_PQR, dir_set, group_set, zonej_set, zonek_set);
  Set const &jplane_set = data_store.newVariable<ProductSet<4>>("Set/JPlane", pspace, SPACE_PQR, dir_set, group_set, zonei_set, zonek_set);
  Set const &kplane_set = data_store.newVariable<ProductSet<4>>("Set/KPlane", pspace, SPACE_PQR, dir_set, group_set, zonei_set, zonej_set);
  data_store.newVariable<Field_IPlane>("i_plane", iplane_set);
  data_store.newVariable<Field_JPlane>("j_plane", jplane_set);
  data_store.newVariable<Field_KPlane>("k_plane", kplane_set);

  // Create a set to span scattering transfer matrix
  Set const &material_set   = data_store.getVariable<Set>("Set/Material");
  Set const &legendre_set   = data_store.getVariable<Set>("Set/Legendre");
  Set const &global_group_set = data_store.getVariable<Set>("Set/GlobalGroup");
  ProductSet<4> *sigs_set = new ProductSet<4>(pspace, SPACE_NULL,
      material_set, legendre_set, global_group_set, global_group_set);

  data_store.addVariable("Set/SigmaS", sigs_set);


  // Create storage for the scattering transfer matrix
  data_store.addVariable("data/sigs", new Field_SigmaS(*sigs_set));
  auto &field_sigs = data_store.getVariable<Field_SigmaS>("data/sigs");

  // Assign basic diagonal data to matrix
  for(auto sdom_id : field_sigs.getWorkList()){

    // Zero out entire matrix
    auto sigs_ptr = field_sigs.getData(sdom_id);
    size_t sigs_size = field_sigs.size(sdom_id);
    for(size_t i = 0;i < sigs_size;++ i){
      sigs_ptr[i] = 0.0;
    }

    // Assign diagonal to the user input for each material
    // Assume each group has same behavior
    auto sigs = field_sigs.getView(sdom_id);
    int global_num_groups = global_group_set.size(sdom_id);
    Legendre n{0};
    for(Material mat{0};mat < 3;++ mat){
      for(GlobalGroup g{0};g < global_num_groups;++ g){
        sigs(mat, n, g, g) = input_vars.sigs[*mat];
      }
    }
  }

  Comm default_comm;
  if(default_comm.rank() == 0){

    unsigned long flux_size = flux_set->globalSize();
    unsigned long moment_size = fluxmoment_set->globalSize();
    unsigned long sigs_size = sigs_set->globalSize();
    unsigned long total_size = (2*moment_size + 2*flux_size + sigs_size);

    printf("\n");
    printf("  Variable     Num Elements      Megabytes\n");
    printf("  --------     ------------      ---------\n");
    printf("  psi          %12lu   %12.3lf\n",
        flux_size,
        (double)flux_size*8.0/1024.0/1024.0);

    printf("  rhs          %12lu   %12.3lf\n",
        flux_size,
        (double)flux_size*8.0/1024.0/1024.0);

    printf("  phi          %12lu   %12.3lf\n",
        moment_size,
        (double)moment_size*8.0/1024.0/1024.0);

    printf("  phi_out      %12lu   %12.3lf\n",
        moment_size,
        (double)moment_size*8.0/1024.0/1024.0);

    printf("  sigs         %12lu   %12.3lf\n",
        sigs_size,
        (double)sigs_size*8.0/1024.0/1024.0);

    printf("  TOTAL        %12lu   %12.3lf\n",
        total_size,
        (double)total_size*8.0/1024.0/1024.0);
  }

}





void Kripke::generateProblem(Kripke::DataStore &data_store,
    InputVariables const &input_vars)
{

  Comm default_comm;

  if(default_comm.rank() == 0){
    printf("\nGenerating\n");
    printf("==========\n\n");
  }

  // Create and start a timing object
  data_store.addVariable("timing", new Kripke::Timing());
  KRIPKE_TIMER(data_store, Generate);


  // Create parallel and subdomain decomposition
  initializeDecomp(data_store, input_vars);

  // Create energy discretization
  initializeEnergy(data_store, input_vars);

  // Create angular discretization, quadrature set and L/L+ matrices
  initializeDirections(data_store, input_vars);

  // Create a spatial mesh, and paint it with materials
  initializeSpace(data_store, input_vars);

  // Create cross sections and transfer matrix
  initializeData(data_store, input_vars);
}
