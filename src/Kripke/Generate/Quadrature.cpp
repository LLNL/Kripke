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

  /**
   * Contains information needed for one quadrature set direction.
   */
  struct QuadraturePoint {
    double xcos;              /* Absolute value of the x-direction cosine. */
    double ycos;              /* Absolute value of the y-direction cosine. */
    double zcos;              /* Absolute value of the z-direction cosine. */
    double w;                 /* weight for the quadrature rule.*/
    int id;                   /* direction flag (= 1 if x-direction
                              cosine is positive; = -1 if not). */
    int jd;                   /* direction flag (= 1 if y-direction
                              cosine is positive; = -1 if not). */
    int kd;                   /* direction flag (= 1 if z-direction
                              cosine is positive; = -1 if not). */
    int octant;
  };


  /*
    GaussLegendre returns the n point Gauss-Legendre quadrature rule for
    the integral between x1 and x2.
  */
  void GaussLegendre(double x1, double x2, std::vector<double> &x,
      std::vector<double> &w, double eps)
  {
    int n = x.size();
    int m, j, i;
    double z1, z, xm, xl, pp, p3, p2, p1;

    m=(n+1)/2;
    xm=0.5*(x2+x1);
    xl=0.5*(x2-x1);
    for(i=1; i<=m; i++){
      z=cos(M_PI*(i-0.25)/(n+0.5));
      do {
        p1=1.0;
        p2=0.0;
        for(j=1; j<=n; j++){
          p3=p2;
          p2=p1;
          p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
        }
        pp=n*(z*p1-p2)/(z*z-1.0);
        z1=z;
        z=z1-p1/pp;
      } while(fabs(z-z1) > eps);
      x[i-1]=xm-xl*z;
      x[n-i]=xm+xl*z;
      w[i-1]=2.0*xl/((1.0-z*z)*pp*pp);

      w[n-i]=w[i-1];
    }
  }


  bool dirSortFcn(QuadraturePoint const &a, QuadraturePoint const &b){
    return b.octant < a.octant;
  }

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
      KRIPKE_ABORT("Bad input to PnmFcn: abs(x) > 1.0, x = %e\n", x);
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



/**
 * Initializes the quadrature set information for a Grid_Data object.
 * This guarantees that each <GS,DS> pair have a single originating octant.
 */
static
std::vector<QuadraturePoint>
createQuadratureSet(InputVariables const &input_vars)
{
  std::vector<QuadraturePoint> directions;

  // Get set description from user
  int num_directions_per_octant = input_vars.num_directions/8;
  int num_directions = input_vars.num_directions;

  // allocate storage
  directions.resize(num_directions);

  // Are we running a REAL quadrature set?
  int num_polar = input_vars.quad_num_polar;
  int num_azimuth = input_vars.quad_num_azimuthal;

  std::vector<double> polar_cos;
  std::vector<double> polar_weight;
  if(num_polar > 0){
    // make sure the user specified the correct number of quadrature points
    KRIPKE_ASSERT(num_polar % 4 == 0,
        "Must have number of polar angles be a multiple of 4\n");

    KRIPKE_ASSERT(num_azimuth % 2 == 0,
      "Must have number of azimuthal angles be a multiple of 2\n");

    KRIPKE_ASSERT(num_polar*num_azimuth == num_directions,
      "You need to specify %d total directions, not %d\n",
          num_polar*num_azimuth, num_directions);

    // Compute gauss legendre weights
    polar_cos.resize(num_polar);
    polar_weight.resize(num_polar);
    GaussLegendre(-1.0, 1.0, polar_cos, polar_weight, DBL_EPSILON);

    // compute azmuhtal angles and weights
    std::vector<double> az_angle(num_azimuth);
    std::vector<double> az_weight(num_azimuth);
    double dangle = 2.0*M_PI/((double) num_azimuth);

    for(int i=0; i<num_azimuth; i++){
      if(i == 0){
        az_angle[0] = dangle/2.0;
      }
      else{
        az_angle[i] = az_angle[i-1] + dangle;
      }
      az_weight[i] = dangle;
    }


    // Loop over polar 'octants
    int d = 0;
    for(int i=0; i< num_polar; i++){
      for(int j=0; j< num_azimuth; j++){
        double xcos = sqrt(1.0-polar_cos[i]*polar_cos[i]) * cos(az_angle[j]);
        double ycos = sqrt(1.0-polar_cos[i]*polar_cos[i]) * sin(az_angle[j]);
        double zcos = polar_cos[i];
        double w = polar_weight[i]*az_weight[j];

        directions[d].id = (xcos > 0.) ? 1 : -1;
        directions[d].jd = (ycos > 0.) ? 1 : -1;
        directions[d].kd = (zcos > 0.) ? 1 : -1;

        directions[d].octant = 0;
        if(directions[d].id == -1){
          directions[d].octant += 1;
        }
        if(directions[d].jd == -1){
          directions[d].octant += 2;
        }
        if(directions[d].kd == -1){
          directions[d].octant += 4;
        }

        directions[d].xcos = std::abs(xcos);
        directions[d].ycos = std::abs(ycos);
        directions[d].zcos = std::abs(zcos);
        directions[d].w = w;

        ++ d;
      }
    }

    // Sort by octant.. so each set has same directions
    std::sort(directions.begin(), directions.end(), dirSortFcn);
  }
  else{
    // Do (essentialy) an S2 quadrature.. but with repeated directions

    // Compute x,y,z cosine values
    double mu  = cos(M_PI/4);
    double eta = sqrt(1-mu*mu) * cos(M_PI/4);
    double xi  = sqrt(1-mu*mu) * sin(M_PI/4);
    int d = 0;
    for(int octant = 0;octant < 8;++ octant){
      double omegas[3];
      omegas[0] = octant & 0x1;
      omegas[1] = (octant>>1) & 0x1;
      omegas[2] = (octant>>2) & 0x1;

      for(int sd=0; sd<num_directions_per_octant; sd++, d++){
        // Store which logical direction of travel we have
        directions[d].id = (omegas[0] > 0.) ? 1 : -1;
        directions[d].jd = (omegas[1] > 0.) ? 1 : -1;
        directions[d].kd = (omegas[2] > 0.) ? 1 : -1;

        // Store quadrature point's weight
        directions[d].w = 4.0*M_PI / (double)num_directions;
        directions[d].xcos = mu;
        directions[d].ycos = eta;
        directions[d].zcos = xi;
      }
    }
  }

  return directions;
}




void Kripke::Generate::generateQuadrature(Kripke::Core::DataStore &data_store,
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
  for(SdomId sdom_id : field_moment_to_legendre.getWorkList()){
    auto moment_to_legendre = field_moment_to_legendre.getView(sdom_id);

    for(Moment nm{0};nm < moment_set->size(sdom_id);++ nm){
      moment_to_legendre(nm) = moment_list[(*nm) + moment_set->lower(sdom_id)];
    }
  }


  // Create the quadrature set
  auto quadrature_points = createQuadratureSet(input_vars);

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
  
	auto &set_ell_plus = data_store.newVariable<ProductSet<2>>("Set/EllPlus",
      pspace, SPACE_Q, *dir_set, *moment_set);

  // Allocate and initialize the L and L+ matrices
  auto &field_ell = data_store.newVariable<Field_Ell>("ell", set_ell);
  auto &field_ell_plus = data_store.newVariable<Field_EllPlus>("ell_plus", set_ell_plus);

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
          ell_plus(d,nm) = ynm*SQRT4PI;
        }
        nm ++;
      }
    }
  }


  // Create fields to store subdomain adjacency information for boundary comm
  // used in sweeps and block jacobi solves
  auto &set_dimension =
      data_store.newVariable<GlobalRangeSet>("Set/Dimension", pspace, 3);

  auto &set_adjacency = data_store.newVariable<ProductSet<1>>("Set/Adjacency",
      pspace, SPACE_PQR, set_dimension);

  auto &field_upwind = data_store.newVariable<Field_Adjacency>("upwind", set_adjacency);
  auto &field_downwind = data_store.newVariable<Field_Adjacency>("downwind", set_adjacency);

  for(SdomId sdom_id : field_upwind.getWorkList()){
    // Get local subdomain coordinates
    auto local_coord = pspace.sdomIdToCoord(sdom_id);

    // Offset local coordinate to global coordinates
    auto global_coord = pspace.coordToGlobalCoord(local_coord);

    std::array<Field_Direction2Int::DefaultViewType, 3> sweep_dir =
        {{
          field_id.getView(sdom_id),
          field_jd.getView(sdom_id),
          field_kd.getView(sdom_id)
        }};

    // Compute upwind and downwind coordinate
    auto upwind = field_upwind.getView(sdom_id);
    auto downwind = field_downwind.getView(sdom_id);
    for(Dimension dim{0};dim < 3;++ dim){

      // Compute upwind and downwind coordinates for this dimensions neighbor
      auto global_upwind = global_coord;
      auto global_downwind = global_coord;
      global_upwind  [*dim+SPACE_RX] -= sweep_dir[*dim](Direction{0});
      global_downwind[*dim+SPACE_RX] += sweep_dir[*dim](Direction{0});

      // Is this an upwind boundary condition?
      if(global_upwind[*dim+SPACE_RX] < 0 ||
         global_upwind[*dim+SPACE_RX] >= (ptrdiff_t)pspace.getGlobalNumSubdomains((Kripke::Core::SPACE)(*dim+SPACE_RX)))
      {
        upwind(dim) = GlobalSdomId{-1};
      }
      // Not a BC, so compute the subdomain id
      else{
        upwind(dim) = pspace.coordToGlobalSdomId(global_upwind);
      }


      // Is this an downwind boundary condition?
      if(global_downwind[*dim+SPACE_RX] < 0 ||
         global_downwind[*dim+SPACE_RX] >= (ptrdiff_t)pspace.getGlobalNumSubdomains((Kripke::Core::SPACE)(*dim+SPACE_RX)))
      {
        downwind(dim) = GlobalSdomId{-1};
      }
      // Not a BC, so compute the subdomain id
      else{
        downwind(dim) = pspace.coordToGlobalSdomId(global_downwind);
      }
    }

  }
}


