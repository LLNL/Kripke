#include "../Param.h"

void dgz_nmd_SweepDD_1d(dgz_nmd_Param &p){
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

  dgz_TVec psi_bo_v(p.num_groups, p.num_directions+1, p.num_zones);
  double ***psi_bo = psi_bo_v.data;

  dgz_TVec psi_l_v(p.num_groups, p.num_directions, p.num_zones+1);
  dgz_TVec psi_r_v(p.num_groups, p.num_directions, p.num_zones+1);
  double ***psi_l = psi_l_v.data;
  double ***psi_r = psi_r_v.data;


  /* Do mu = -1 first */
  double **psi_d0 = p.psi.data[0];
  double **RHS_d0 = p.rhs.data[0];
  double **psi_bo_d0 = psi_bo[0];
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(int group = 0; group < p.num_groups; ++group){
    double * __restrict__ psi_d0_g = psi_d0[group];
    double * __restrict__ RHS_d0_g = RHS_d0[group];
    double * __restrict__ psi_bo_d0_g = psi_bo_d0[group];
    double * __restrict__ sigt_g = p.sigt.data[group];

    // pull out first zone
    double psi_d0_g_z0 = psi_d0_g[0];
    double psi_r1 = psi_d0_g_z0;   /* psi[0][0] has isotropic boundary value.
                                     */
    for(int i=num_zones-1; i>=0; i--){
      psi_r1 += psi_internal_in_g[i];   /* Add internal source if any */
      double dxip1 = delta_x[i+1];
      double dxip1_RHS_d0_g_z = dxip1*RHS_d0_g[i];
      double dxip1_sigt_g_z = dxip1*sigt_g[i];
      double psi_bo_0_i =
        (2.0*psi_r1 + dxip1_RHS_d0_g_z)/(2.0 + dxip1_sigt_g_z);
      psi_r1 = 2.0*psi_bo_0_i - psi_r1;
      psi_bo_d0_g[i] = psi_bo_0_i;
    }
    for(int d=0; d<Lm; d++){
      /* Load values at x=0. */
      /*psi_r[d][0] = psi_r1;*/
      /* psi[0][0] has isotropic boundary value. */
      psi_r[d][group][num_zones] = psi_d0_g_z0;
    }
    for(int d=0; d<num_directions; d++){
      /* Load values at x=0. */
      psi_l[d][group][0] = psi_r1;

    }
  }

  /* Do mu < 0 next */
  for(int d=0; d<Lm; d++){
    double **psi_d = p.psi.data[d];
    double **RHS_d = p.rhs.data[d];
    double **psi_bo_d = psi_bo[d];
    double **psi_bo_dp1 = psi_bo[d+1];
    double **psi_r_d = psi_r[d];
    double **psi_l_d = psi_l[d];

    double aup = direction[d].alpha_p;
    double adn = direction[d].alpha_m;
    double mu = direction[d].xcos;
    double w = direction[d].w;

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(int group = 0; group < p.num_groups; ++group){
      double * __restrict__ psi_d_g = psi_d[group];
      double * __restrict__ RHS_d_g = RHS_d[group];
      double * __restrict__ psi_bo_d_g = psi_bo_d[group];
      double * __restrict__ psi_bo_dp1_g = psi_bo_dp1[group];
      double * __restrict__ psi_r_d_g = psi_r_d[group];
      double * __restrict__ psi_l_d_g = psi_l_d[group];
      double * __restrict__ sigt_g = p.sigt.data[group];

      for(int i=num_zones-1; i>=1; i--){
        double Ai = area[i];
        double Aip1 = area[i+1];
        double Vi = volume[i];

        psi_r_d_g[i+1] += psi_internal_in_g[i]; /* Add internal source if any
                                                  */

        double A1_tmp = mu*(Aip1 + Ai);
        double A2_tmp = ((Aip1 - Ai)/w)*(aup+adn);
        double A3_tmp = 2.0*(Aip1 - Ai)*aup/w;
        double psi_d_i = 2.0 *
                         ( A1_tmp*psi_r_d_g[i+1] + A2_tmp*psi_bo_d_g[i]
                           + Vi*
                           RHS_d_g[i] ) / ( Vi*sigt_g[i] + 2.0*mu*Ai + A3_tmp );
        psi_r_d_g[i] = psi_d_i - psi_r_d_g[i+1];
        psi_bo_dp1_g[i] = psi_d_i - psi_bo_d_g[i];
        psi_d_g[i] = 0.5*psi_d_i;
      }

      /* do i=0 */
      double A0 = area[0];
      double A1 = area[1];
      double V0 = volume[0];
      double psi_d_g_0 =
        ( mu*(A1*psi_r_d_g[1] - A0*psi_r_d_g[0]) +
          (A1-A0)/(w)*(aup+adn)*psi_bo_d_g[0]
          + V0*RHS_d_g[0] ) / ( V0*sigt_g[0] + 2.0*(A1-A0)*aup/w );

      psi_bo_dp1_g[0] = 2.0*psi_d_g_0 - psi_bo_d_g[0];
      psi_d_g[0] = psi_d_g_0;
    }
  }

  /* Do mu > 0 last */
  for(int d=Lm; d<num_directions; d++){
    double **psi_d = p.psi.data[d];
    double **RHS_d = p.rhs.data[d];
    double **psi_bo_d = psi_bo[d];
    double **psi_bo_dp1 = psi_bo[d+1];
    double **psi_r_d = psi_r[d];
    double **psi_l_d = psi_l[d];

    double aup = direction[d].alpha_p;
    double adn = direction[d].alpha_m;
    double mu = direction[d].xcos;
    double w = direction[d].w;

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(int group = 0; group < p.num_groups; ++group){
      double * __restrict__ psi_d_g = psi_d[group];
      double * __restrict__ RHS_d_g = RHS_d[group];
      double * __restrict__ psi_bo_d_g = psi_bo_d[group];
      double * __restrict__ psi_bo_dp1_g = psi_bo_dp1[group];
      double * __restrict__ psi_r_d_g = psi_r_d[group];
      double * __restrict__ psi_l_d_g = psi_l_d[group];
      double * __restrict__ sigt_g = p.sigt.data[group];
      for(int i=0; i<num_zones; i++){
        double Ai = area[i];
        double Aip1 = area[i+1];
        double Vi = volume[i];

        psi_l_d_g[i] += psi_internal_out_g[i]; /* Add internal source if any
                                                 */

        double A1_tmp = mu*(Aip1 + Ai);
        double A2_tmp = ((Aip1 - Ai)/w)*(aup+adn);
        double A3_tmp = 2.0*(Aip1 - Ai)*aup/w;
        double psi_d_g_i = 2.0 *
                           ( A1_tmp*psi_l_d_g[i] + A2_tmp*psi_bo_d_g[i]
                             + Vi*
                             RHS_d_g[i] ) /
                           ( Vi*sigt_g[i] + 2.0*mu*Aip1 + A3_tmp );
        psi_l_d_g[i+1] = psi_d_g_i - psi_l_d_g[i];
        psi_bo_dp1_g[i] = psi_d_g_i - psi_bo_d_g[i];
        psi_d_g[i] = 0.5*psi_d_g_i;
      }
    }
  }

  /* Load psi_out */
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(int group = 0; group < p.num_groups; ++group){
    for(int d=0; d<Lm; d++){
      psi_out[group][d] = 0.0;
    }
    for(int d=Lm; d<num_directions; d++){
      psi_out[group][d] = psi_l[d][group][num_zones];
    }
  }

}



void dgz_nmd_SweepDD_3d(dgz_nmd_Param &p);

void dgz_nmd_Param::SweepDD(void) {

  if(geometry_type == 1){
    dgz_nmd_SweepDD_1d(*this);
  }
  else if(geometry_type == 2){
    //dgz_nmd_SweepDD_2d(*this);
  }
  else {
    dgz_nmd_SweepDD_3d(*this);
  }
}
