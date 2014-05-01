#include "../Param.h"
#include <vector>

void dzg_nmd_SweepDD_1d(dzg_nmd_Param &p){
  int num_directions = p.num_directions;
  int num_zones = p.num_zones;
  Direction *direction = p.direction;
  double *volume = p.volume;
  double *area = p.area;
  double *delta_x = p.deltas[0];  /* Ghosted array */

  // do these need to be per-group?
  double *psi_internal_in_g = p.tmp_source_in;
  double *psi_internal_out_g= p.tmp_source_out;

  /* Calculate Lp=no. of positive mu values, and Lm = num_directions - Lp */
  int Lp = 0;
  for(int d=0; d<num_directions; d++){
    if(direction[d].id > 0){
      Lp++;
    }
  }
  int Lm = num_directions - Lp;

  /* Allocate temporary storage */
  gz_Vec psi_out_v(p.num_groups, p.num_directions);
  double **psi_out = psi_out_v.data;

  dzg_TVec psi_bo_v(p.num_groups, p.num_directions+1, p.num_zones);
  double ***psi_bo = psi_bo_v.data;

  dzg_TVec psi_l_v(p.num_groups, p.num_directions, p.num_zones+1);
  dzg_TVec psi_r_v(p.num_groups, p.num_directions, p.num_zones+1);
  double ***psi_l = psi_l_v.data;
  double ***psi_r = psi_r_v.data;

  std::vector<double> psi_r1_v;
  psi_r1_v.resize(p.num_groups);
  double * __restrict__ psi_r1 = &(psi_r1_v[0]);

  std::vector<double> psi_d0_z0_v;
  psi_d0_z0_v.resize(p.num_groups);
  double * __restrict__ psi_d0_z0 = &(psi_d0_z0_v[0]);

  /* Do mu = -1 first */
  double **psi_d0 = p.psi.data[0];
  double **RHS_d0 = p.rhs.data[0];
  double **psi_bo_d0 = psi_bo[0];
  for(int group = 0; group < p.num_groups; ++group){
    psi_r1[group] = psi_d0[0][group];
    psi_d0_z0[group] = psi_d0[0][group];
  }
  for(int i=num_zones-1; i>=0; i--){
    double * __restrict__ RHS_d0_z = RHS_d0[i];
    double * __restrict__ psi_bo_d0_z = psi_bo_d0[i];
    double * __restrict__ sigt_z = p.sigt.data[i];

    double dxip1 = delta_x[i+1];

    for(int group = 0; group < p.num_groups; ++group){
      // pull out first zone
      double psi_d0_z0_g = psi_d0_z0[group];

      psi_r1[group] += psi_internal_in_g[i];   /* Add internal source if any
                                                 */
      double dxip1_RHS_d0_g_z = dxip1*RHS_d0_z[group];
      double dxip1_sigt_g_z = dxip1*sigt_z[group];
      double psi_bo_d0_z_g =
        (2.0*psi_r1[group] + dxip1_RHS_d0_g_z)/(2.0 + dxip1_sigt_g_z);
      psi_r1[group] = 2.0*psi_bo_d0_z_g - psi_r1[group];
      psi_bo_d0_z[group] = psi_bo_d0_z_g;

    }
  }
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(int group = 0; group < p.num_groups; ++group){
    for(int d=0; d<Lm; d++){
      /* Load values at x=0. */
      /*psi_r[d][0] = psi_r1;*/
      /* psi[0][0] has isotropic boundary value. */
      psi_r[d][num_zones][group] = psi_d0_z0[group];
    }
    for(int d=0; d<num_directions; d++){
      /* Load values at x=0. */
      psi_l[d][0][group] = psi_r1[group];
    }
  }
  /* Do mu < 0 next */
  for(int d=0; d<Lm; d++){
    double **psi_d = p.psi.data[d];
    double **RHS_d = p.rhs.data[d];
    double **psi_bo_d = psi_bo[d];
    double **psi_bo_dp1 = psi_bo[d+1];
    double **psi_r_d = psi_r[d];

    double aup = direction[d].alpha_p;
    double adn = direction[d].alpha_m;
    double mu = direction[d].xcos;
    double w = direction[d].w;

    for(int i=num_zones-1; i>=1; i--){
      double Ai = area[i];
      double Aip1 = area[i+1];
      double Vi = volume[i];
      double A1_tmp = mu*(Aip1 + Ai);
      double A2_tmp = ((Aip1 - Ai)/w)*(aup+adn);
      double A3_tmp = 2.0*(Aip1 - Ai)*aup/w;

      double * __restrict__ psi_d_z = psi_d[i];
      double * __restrict__ RHS_d_z = RHS_d[i];
      double * __restrict__ psi_bo_d_z = psi_bo_d[i];
      double * __restrict__ psi_bo_dp1_z = psi_bo_dp1[i];
      double * __restrict__ psi_r_d_z = psi_r_d[i];
      double * __restrict__ psi_r_d_zp1 = psi_r_d[i+1];
      double * __restrict__ sigt_z = p.sigt.data[i];

      for(int group = 0; group < p.num_groups; ++group){
        psi_r_d_zp1[group] += psi_internal_in_g[i]; /* Add internal source if any
                                                  */

        double psi_d_z_g = 2.0 *
                           ( A1_tmp*psi_r_d_zp1[group] + A2_tmp*
                             psi_bo_d_z[group]
                             + Vi*
                             RHS_d_z[group] ) /
                           ( Vi*sigt_z[group] + 2.0*mu*Ai + A3_tmp );
        psi_r_d_z[group] = psi_d_z_g - psi_r_d_zp1[group];
        psi_bo_dp1_z[group] = psi_d_z_g - psi_bo_d_z[group];
        psi_d_z[group] = 0.5*psi_d_z_g;
      }
    }

    /* do i=0 */
    double * __restrict__ psi_r_d_z0 = psi_r_d[0];
    double * __restrict__ psi_d_z0 = psi_d[0];
    double * __restrict__ psi_r_d_z1 = psi_r_d[1];
    double * __restrict__ psi_bo_d_z0 = psi_bo_d[0];
    double * __restrict__ psi_bo_dp1_z0 = psi_bo_dp1[0];
    double * __restrict__ RHS_d_z0 = RHS_d[0];
    double * __restrict__ sigt_z0 = p.sigt.data[0];
    for(int group = 0; group < p.num_groups; ++group){
      double A0 = area[0];
      double A1 = area[1];
      double V0 = volume[0];
      double psi_d_z0_g =
        ( mu*(A1*psi_r_d_z1[group] - A0*psi_r_d_z0[group]) +
          (A1-A0)/(w)*(aup+adn)*psi_bo_d_z0[group]
          + V0*RHS_d_z0[group] ) / ( V0*sigt_z0[group] + 2.0*(A1-A0)*aup/w );

      psi_bo_dp1_z0[group] = 2.0*psi_d_z0_g - psi_bo_d_z0[group];
      psi_d_z0[group] = psi_d_z0_g;
    }
  }


  /* Do mu > 0 last */
  for(int d=Lm; d<num_directions; d++){
    double **psi_d = p.psi.data[d];
    double **RHS_d = p.rhs.data[d];
    double **psi_bo_d = psi_bo[d];
    double **psi_bo_dp1 = psi_bo[d+1];
    double **psi_l_d = psi_l[d];

    double aup = direction[d].alpha_p;
    double adn = direction[d].alpha_m;
    double mu = direction[d].xcos;
    double w = direction[d].w;

    for(int i=0; i<num_zones; i++){
      double Ai = area[i];
      double Aip1 = area[i+1];
      double Vi = volume[i];
      double A1_tmp = mu*(Aip1 + Ai);
      double A2_tmp = ((Aip1 - Ai)/w)*(aup+adn);
      double A3_tmp = 2.0*(Aip1 - Ai)*aup/w;

      double * __restrict__ psi_d_z = psi_d[i];
      double * __restrict__ RHS_d_z = RHS_d[i];
      double * __restrict__ psi_bo_d_z = psi_bo_d[i];
      double * __restrict__ psi_bo_dp1_z = psi_bo_dp1[i];
      double * __restrict__ psi_l_d_z = psi_l_d[i];
      double * __restrict__ psi_l_d_zp1 = psi_l_d[i+1];
      double * __restrict__ sigt_z = p.sigt.data[i];

      for(int group = 0; group < p.num_groups; ++group){
        psi_l_d_z[group] += psi_internal_out_g[i]; /* Add internal source if any
                                                  */

        double psi_d_z_g = 2.0 *
                           ( A1_tmp*psi_l_d_z[group] + A2_tmp*
                             psi_bo_d_z[group]
                             + Vi*
                             RHS_d_z[group] ) /
                           ( Vi*sigt_z[group] + 2.0*mu*Aip1 + A3_tmp );
        psi_l_d_zp1[group] = psi_d_z_g - psi_l_d_z[group];
        psi_bo_dp1_z[group] = psi_d_z_g - psi_bo_d_z[group];
        psi_d_z[group] = 0.5*psi_d_z_g;
      }
    }
  }

  /* Load psi_out */
  for(int group = 0; group < p.num_groups; ++group){
    for(int d=0; d<Lm; d++){
      psi_out[group][d] = 0.0;
    }
    for(int d=Lm; d<num_directions; d++){
      psi_out[group][d] = psi_l[d][group][num_zones];
    }
  }
}



void dzg_nmd_SweepDD_3d(dzg_nmd_Param &p);

void dzg_nmd_Param::SweepDD(void) {

  if(geometry_type == 1){
    dzg_nmd_SweepDD_1d(*this);
  }
  else if(geometry_type == 2){
    //dzg_nmd_SweepDD_2d(*this);
  }
  else {
    dzg_nmd_SweepDD_3d(*this);
  }
}
