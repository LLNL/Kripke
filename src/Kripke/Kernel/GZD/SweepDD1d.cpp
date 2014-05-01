#include "../Param.h"

void gzd_nmd_SweepDD_1d(gzd_nmd_Param &p){
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

  gzd_TVec psi_bo_v(p.num_groups, p.num_directions+1, p.num_zones);
  double ***psi_bo = psi_bo_v.data;

  gzd_TVec psi_l_v(p.num_groups, p.num_directions, p.num_zones+1);
  gzd_TVec psi_r_v(p.num_groups, p.num_directions, p.num_zones+1);
  double ***psi_l = psi_l_v.data;
  double ***psi_r = psi_r_v.data;


#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(int group = 0; group < p.num_groups; ++group){

    double **psi_bo_g = psi_bo[group];
    double **psi_l_g = psi_l[group];
    double **psi_r_g = psi_r[group];
    double * __restrict__ psi_out_g = psi_out[group];

    double **psi_g = p.psi.data[group];
    double **RHS_g = p.rhs.data[group];
    double *sigt_g = p.sigt.data[group];
    for(int z=0; z<num_zones; z++){
      std::fill_n(psi_bo_g[z], num_directions+1, 0.0);
    }
    for(int z=0; z<num_zones+1; z++){
      std::fill_n(psi_l_g[z], num_directions, 0.0);
      std::fill_n(psi_r_g[z], num_directions, 0.0);
    }


    /* Do mu = -1 first */
    double psi_0_0 = psi_g[0][0];
    double psi_r_g1 = psi_0_0;   /* psi[0][0] has isotropic boundary value. */
    for(int i=num_zones-1; i>=0; i--){
      psi_r_g1 += psi_internal_in_g[i];   /* Add internal source if any */
      double dxip1 = delta_x[i+1];
      double dxip1_RHS_i_0 = dxip1*RHS_g[i][0];
      double dxip1_sigt_g_z = dxip1*sigt_g[i];
      double psi_bo_g_0_i =
        (2.0*psi_r_g1 + dxip1_RHS_i_0)/(2.0 + dxip1_sigt_g_z);
      psi_r_g1 = 2.0*psi_bo_g_0_i - psi_r_g1;
      psi_bo_g[i][0] = psi_bo_g_0_i;
    }
    for(int d=0; d<Lm; d++){
      /* Load values at x=0. */
      /*psi_r_g[d][0] = psi_r_g1;*/
      /* psi[0][0] has isotropic boundary value. */
      psi_r_g[num_zones][d] = psi_0_0;
    }
    for(int d=0; d<num_directions; d++){
      /* Load values at x=0. */
      psi_l_g[0][d] = psi_r_g1;
    }


    /* Do mu < 0 next */
    for(int i=num_zones-1; i>=1; i--){
      double Ai = area[i];
      double Aip1 = area[i+1];
      double Vi = volume[i];

      double * __restrict__ psi_g_z = psi_g[i];
      double * __restrict__ psi_r_g_z = psi_r_g[i];
      double * __restrict__ psi_r_g_zp1 = psi_r_g[i+1];
      double * __restrict__ psi_bo_g_z = psi_bo_g[i];
      double * __restrict__ RHS_g_z = RHS_g[i];

      for(int d=0; d<Lm; d++){

        double aup = direction[d].alpha_p;
        double adn = direction[d].alpha_m;
        double mu = direction[d].xcos;
        double w = direction[d].w;
        //double winv = 1.0/w;

        psi_r_g_zp1[d] += psi_internal_in_g[i]; /* Add internal source if any
                                                  */

        double A1_tmp = mu*(Aip1 + Ai);
        double A2_tmp = ((Aip1 - Ai)/w)*(aup+adn);
        double A3_tmp = 2.0*(Aip1 - Ai)*aup/w;
        double psi_d_i = 2.0 *
                         ( A1_tmp*psi_r_g_zp1[d] + A2_tmp*psi_bo_g_z[d]
                           + Vi*
                           RHS_g_z[d] ) / ( Vi*sigt_g[i] + 2.0*mu*Ai + A3_tmp );
        psi_r_g_z[d] = psi_d_i - psi_r_g_zp1[d];
        psi_bo_g_z[d+1] = psi_d_i - psi_bo_g_z[d];
        psi_g_z[d] = 0.5*psi_d_i;
      }
    }

    /* do i=0 */
    {
      double * __restrict__ psi_r_g_z0 = psi_r_g[0];
      double * __restrict__ psi_r_g_z1 = psi_r_g[1];

      double * __restrict__ psi_bo_g_z0 = psi_bo_g[0];
      double * __restrict__ F_z0 = RHS_g[0];
      double * __restrict__ psi_g_z0 = psi_g[0];
      double A0 = area[0];
      double A1 = area[1];
      double V0 = volume[0];
      for(int d=0; d<Lm; d++){
        double aup = direction[d].alpha_p;
        double adn = direction[d].alpha_m;
        double mu = direction[d].xcos;
        double w = direction[d].w;

        double psi_d_0 =
          ( mu*(A1*psi_r_g_z1[d] - A0*psi_r_g_z0[d]) +
            (A1-A0)/(w)*(aup+adn)*psi_bo_g_z0[d]
            + V0*F_z0[d] ) / ( V0*sigt_g[0] + 2.0*(A1-A0)*aup/w );

        psi_bo_g_z0[d+1] = 2.0*psi_d_0 - psi_bo_g_z0[d];
        psi_g_z0[d] = psi_d_0;
      }
    }

    /* Do mu > 0 last */
    for(int i=0; i<num_zones; i++){
      double Ai = area[i];
      double Aip1 = area[i+1];
      double Vi = volume[i];

      double * __restrict__ psi_g_z = psi_g[i];
      double * __restrict__ psi_l_g_z = psi_l_g[i];
      double * __restrict__ psi_l_g_zp1 = psi_l_g[i+1];
      double * __restrict__ psi_bo_g_z = psi_bo_g[i];
      double * __restrict__ RHS_g_z = RHS_g[i];

      for(int d=Lm; d<num_directions; d++){

        double aup = direction[d].alpha_p;
        double adn = direction[d].alpha_m;
        double mu = direction[d].xcos;
        double w = direction[d].w;

        psi_l_g_z[d] += psi_internal_out_g[i]; /* Add internal source if any
                                                 */

        double A1_tmp = mu*(Aip1 + Ai);
        double A2_tmp = ((Aip1 - Ai)/w)*(aup+adn);
        double A3_tmp = 2.0*(Aip1 - Ai)*aup/w;
        double psi_d_i = 2.0 *
                         ( A1_tmp*psi_l_g_z[d] + A2_tmp*psi_bo_g_z[d]
                           + Vi*
                           RHS_g_z[d] ) /
                         ( Vi*sigt_g[i] + 2.0*mu*Aip1 + A3_tmp );
        psi_l_g_zp1[d] = psi_d_i - psi_l_g_z[d];
        psi_bo_g_z[d+1] = psi_d_i - psi_bo_g_z[d];
        psi_g_z[d] = 0.5*psi_d_i;
      }
    }

    /* Load psi_out_g */
    for(int d=0; d<Lm; d++){
      psi_out_g[d] = 0.0;
    }
    for(int d=Lm; d<num_directions; d++){
      psi_out_g[d] = psi_l_g[num_zones][d];
    }


  }


}


void gzd_nmd_SweepDD_3d(gzd_nmd_Param &p);

void gzd_nmd_Param::SweepDD(void) {

  if(geometry_type == 1){
    gzd_nmd_SweepDD_1d(*this);
  }
  else if(geometry_type == 2){
    //gzd_nmd_SweepDD_2d(*this);
  }
  else {
    gzd_nmd_SweepDD_3d(*this);
  }
}
