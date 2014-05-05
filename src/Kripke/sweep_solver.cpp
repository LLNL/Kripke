/*--------------------------------------------------------------------------
 * Sweep-based solver routine.
 *--------------------------------------------------------------------------*/

#include<Kripke/comm.h>
#include<Kripke/user_data.h>
#include<Kripke/transport_protos.h>
#include<vector>
#include<stdio.h>


/* Local prototypes */
void CreateBufferInfoDD3D(User_Data *user_data);
int SweepSolverSolveDD (User_Data *user_data, double **rhs,
                        double **ans, double **tempv);


/*----------------------------------------------------------------------
 * SweepSolverSolve
 *----------------------------------------------------------------------*/

int SweepSolverSolve (User_Data *user_data, double **rhs, double **ans,
                      double **tempv)
{
  /* Begin timing */
  user_data->timing.start("Sweep");

  /* Diamond Difference */
  SweepSolverSolveDD(user_data, rhs, ans, tempv);

  /* End timing */
  user_data->timing.stop("Sweep");

  return(0);
}


/*----------------------------------------------------------------------
 * SweepSolverSolveDD
 *----------------------------------------------------------------------*/

int SweepSolverSolveDD (User_Data *user_data, double **rhs,
                        double **ans, double **tempv)
{
  Grid_Data  *grid_data         = user_data->grid_data;
  std::vector<Directions> &directions        = user_data->directions;
  std::vector<double>      &tmp_sigma_tot     = grid_data->tmp_sigma_tot;

  int *nzones          = grid_data->nzones;
  int num_zones = grid_data->num_zones;
  int nx=nzones[0], ny=nzones[1], nz=nzones[2];
  int num_directions = user_data->directions.size();

  double *msg;
  int i, k, sweep_group, i_plane_zones, j_plane_zones, k_plane_zones,
      i_src_subd, j_src_subd, k_src_subd, i_dst_subd,
      j_dst_subd, k_dst_subd, directions_left, in, ip, jn, jp, kn, kp, d;
  int local_imax, local_jmax, local_kmax;

  int bc_ref_in, bc_ref_ip, bc_ref_jn, bc_ref_jp, bc_ref_kn, bc_ref_kp;
  double eta_ref_in, eta_ref_ip, eta_ref_jn, eta_ref_jp;
  double eta_ref_kn, eta_ref_kp;
  double eta_ref_i, eta_ref_j, eta_ref_k;

  /*spectral reflection rules relating eminating and iminating fluxes
    for each of the 8 octant for the planes: i,j,k*/
  int r_rules[8][3] = {{1, 3, 4},
                       {0, 2, 5},
                       {3, 1, 6},
                       {2, 0, 7},
                       {5, 7, 0},
                       {4, 6, 1},
                       {7, 5, 2},
                       {6, 4, 3}, };
  int ref_d, octant, ref_octant, fundamental_d;
  int bc_ref_i, bc_ref_j, bc_ref_k;
  int emminating_directions_left;

  bc_ref_in = user_data->bc_types[0];
  bc_ref_ip = user_data->bc_types[1];
  bc_ref_jn = user_data->bc_types[2];
  bc_ref_jp = user_data->bc_types[3];
  bc_ref_kn = user_data->bc_types[4];
  bc_ref_kp = user_data->bc_types[5];

  in = grid_data->mynbr[0][0];
  ip = grid_data->mynbr[0][1];
  jn = grid_data->mynbr[1][0];
  jp = grid_data->mynbr[1][1];
  kn = grid_data->mynbr[2][0];
  kp = grid_data->mynbr[2][1];

  if(in == -1){
    eta_ref_in = user_data->bc_values[0];
  }
  else {
    eta_ref_in = 0.0;
  }
  if(ip == -1){
    eta_ref_ip = user_data->bc_values[1];
  }
  else {
    eta_ref_ip = 0.0;
  }
  if(jn == -1){
    eta_ref_jn = user_data->bc_values[2];
  }
  else {
    eta_ref_jn = 0.0;
  }
  if(jp == -1){
    eta_ref_jp = user_data->bc_values[3];
  }
  else {
    eta_ref_jp = 0.0;
  }
  if(kn == -1){
    eta_ref_kn = user_data->bc_values[4];
  }
  else {
    eta_ref_kn = 0.0;
  }
  if(kp == -1){
    eta_ref_kp = user_data->bc_values[5];
  }
  else {
    eta_ref_kp = 0.0;
  }

  local_imax = nzones[0];
  local_jmax = nzones[1];
  local_kmax = nzones[2];
  i_plane_zones = local_jmax * local_kmax;
  j_plane_zones = local_imax * local_kmax;
  k_plane_zones = local_imax * local_jmax;

  std::vector<double*> i_plane_data(num_directions, NULL);
  std::vector<double*> j_plane_data(num_directions, NULL);
  std::vector<double*> k_plane_data(num_directions, NULL);

  std::vector<int> i_which(num_directions);
  std::vector<int> j_which(num_directions);
  std::vector<int> k_which(num_directions);
  for(d=0; d<num_directions; d++){
    i_which[d] = (directions[d].id>0) ? 0 : 1;
    j_which[d] = (directions[d].jd>0) ? 2 : 3;
    k_which[d] = (directions[d].kd>0) ? 4 : 5;
  }

  std::vector<double> psi_lf_data((local_imax+1)*local_jmax*local_kmax);
  std::vector<double> psi_fr_data( local_imax*(local_jmax+1)*local_kmax);
  std::vector<double> psi_bo_data(local_imax*local_jmax*(local_kmax+1));

  /* Evaluate the total cross section for this group */
  EvalSigmaTot(user_data, tmp_sigma_tot);

  /* Hang out receive requests for each of the 6 neighbors */
  if(in != -1){
    R_recv_dir( 0, in );
  }
  if(ip != -1){
    R_recv_dir( 1, ip );
  }
  if(jn != -1){
    R_recv_dir( 2, jn );
  }
  if(jp != -1){
    R_recv_dir( 3, jp );
  }
  if(kn != -1){
    R_recv_dir( 4, kn );
  }
  if(kp != -1){
    R_recv_dir( 5, kp );
  }

  /* Allocate and initialize (set to zero for now) message
     buffers for subdomain faces on the problem boundary */
  for(d=0; d<num_directions; d++){

    i_src_subd = directions[d].i_src_subd;
    j_src_subd = directions[d].j_src_subd;
    k_src_subd = directions[d].k_src_subd;

    /* get reflective b.c. information for src faces */
    bc_ref_i = (directions[d].id>0) ? bc_ref_in : bc_ref_ip;
    bc_ref_j = (directions[d].jd>0) ? bc_ref_jn : bc_ref_jp;
    bc_ref_k = (directions[d].kd>0) ? bc_ref_kn : bc_ref_kp;
    eta_ref_i = (directions[d].id>0) ? eta_ref_in : eta_ref_ip;
    eta_ref_j = (directions[d].jd>0) ? eta_ref_jn : eta_ref_jp;
    eta_ref_k = (directions[d].kd>0) ? eta_ref_kn : eta_ref_kp;

    if(k_src_subd == -1 && bc_ref_k == 0){
      if(R_recv_test( k_which[d], &(k_plane_data[d]) ) == 0){
        printf("Null buffer not returned to DD_Sweep\n");
        error_exit(1);
      }
      for(k=0; k<k_plane_zones; k++){
        k_plane_data[d][k] = eta_ref_k;
      }
      k_plane_data[d][k_plane_zones] = (double) d;
    }
    else {
      k_plane_data[d] = NULL;
    }

    if(j_src_subd == -1 && bc_ref_j == 0){
      if(R_recv_test( j_which[d], &(j_plane_data[d]) ) == 0){
        printf("Null buffer not returned to DD_Sweep\n");
        error_exit(1);
      }
      for(k=0; k<j_plane_zones; k++){
        j_plane_data[d][k] = eta_ref_j;
      }
      j_plane_data[d][j_plane_zones] = (double) d;
    }
    else {
      j_plane_data[d] = NULL;
    }

    if(i_src_subd == -1 && bc_ref_i == 0){
      if(R_recv_test( i_which[d], &(i_plane_data[d]) ) == 0){
        printf("Null buffer not returned to DD_Sweep\n");
        error_exit(1);
      }
      for(i=0; i<i_plane_zones; i++){
        i_plane_data[d][i] = eta_ref_i;
      }
      i_plane_data[d][i_plane_zones] = (double) d;
    }
    else {
      i_plane_data[d] = NULL;
    }
  }

  directions_left = num_directions;
  std::vector<int> swept(num_directions, 0.0);

  while(directions_left){

    /* Check for a message from the 6 neighboring subdomains. */
    if(in != -1 && R_recv_test( 0, &msg ) != 0){
      i_plane_data[(int) msg[i_plane_zones]] = msg;
    }

    if(ip != -1 && R_recv_test( 1, &msg ) != 0){
      i_plane_data[(int) msg[i_plane_zones]] = msg;
    }

    if(jn != -1 && R_recv_test( 2, &msg ) != 0){
      j_plane_data[(int) msg[j_plane_zones]] = msg;
    }

    if(jp != -1 && R_recv_test( 3, &msg ) != 0){
      j_plane_data[(int) msg[j_plane_zones]] = msg;
    }

    if(kn != -1 && R_recv_test( 4, &msg ) != 0){
      k_plane_data[(int) msg[k_plane_zones]] = msg;
    }

    if(kp != -1 && R_recv_test( 5, &msg ) != 0){
      k_plane_data[(int) msg[k_plane_zones]] = msg;
    }

    for(d=0; d<num_directions; d++){
      if(k_plane_data[d] == NULL ||
         j_plane_data[d] == NULL ||
         i_plane_data[d] == NULL ||
         swept[d]){
        continue;
      }

      /* Use standard Diamond-Difference sweep */
      SweepDD(d, grid_data, user_data->directions, grid_data->volume, tmp_sigma_tot,
              rhs[d], ans[d], i_plane_data[d], j_plane_data[d],
              k_plane_data[d], &psi_lf_data[0], &psi_fr_data[0],
              &psi_bo_data[0]);

      i_dst_subd = directions[d].i_dst_subd;
      j_dst_subd = directions[d].j_dst_subd;
      k_dst_subd = directions[d].k_dst_subd;

      R_send( i_plane_data[d], i_dst_subd, i_plane_zones+1 );
      R_send( j_plane_data[d], j_dst_subd, j_plane_zones+1 );
      R_send( k_plane_data[d], k_dst_subd, k_plane_zones+1 );

      /*Check if any of the 3 planes are reflective problem boundaries.
    If so, generate the src for future sweeps */

      /* get reflective b.c. information for dst faces */
      bc_ref_i = (directions[d].id>0) ? bc_ref_ip : bc_ref_in;
      bc_ref_j = (directions[d].jd>0) ? bc_ref_jp : bc_ref_jn;
      bc_ref_k = (directions[d].kd>0) ? bc_ref_kp : bc_ref_kn;
      eta_ref_i = (directions[d].id>0) ? eta_ref_ip : eta_ref_in;
      eta_ref_j = (directions[d].jd>0) ? eta_ref_jp : eta_ref_jn;
      eta_ref_k = (directions[d].kd>0) ? eta_ref_kp : eta_ref_kn;

      if(k_dst_subd == -1 && bc_ref_k == 1){
        octant = user_data->octant_map[d];
        ref_octant = r_rules[octant][2];
        fundamental_d = (d - octant)/8;
        ref_d = 8 * fundamental_d + ref_octant;
        /* printf("k: d= %2d o=%2d fund_d=%2d ref_o=%2d ref_d=%2d\n",
           d, octant, fundamental_d,ref_octant,ref_d) ; */
        if(R_recv_test( k_which[ref_d], &(k_plane_data[ref_d]) )
           == 0){
          printf("Null buffer not returned to DD_Sweep\n");
          error_exit(1);
        }
        for(k=0; k<k_plane_zones; k++){
          k_plane_data[ref_d][k] = eta_ref_k * k_plane_data[d][k];
        }
        k_plane_data[ref_d][k_plane_zones] = (double) ref_d;
      }
      if(j_dst_subd == -1 && bc_ref_j == 1){
        octant = user_data->octant_map[d];
        ref_octant = r_rules[octant][1];
        fundamental_d = (d - octant)/8;
        ref_d = 8 * fundamental_d + ref_octant;
        /* printf("j: d= %2d o=%2d fund_d=%2d ref_o=%2d ref_d=%2d\n",
           d, octant, fundamental_d,ref_octant,ref_d) ; */
        if(R_recv_test( j_which[ref_d], &(j_plane_data[ref_d]) )
           == 0){
          printf("Null buffer not returned to DD_Sweep\n");
          error_exit(1);
        }
        for(k=0; k<j_plane_zones; k++){
          j_plane_data[ref_d][k] = eta_ref_j * j_plane_data[d][k];
        }
        j_plane_data[ref_d][j_plane_zones] = (double) ref_d;

      }
      if(i_dst_subd == -1 && bc_ref_i == 1){
        octant = user_data->octant_map[d];
        ref_octant = r_rules[octant][0];
        fundamental_d = (d - octant)/8;
        ref_d = 8 * fundamental_d + ref_octant;
        /* printf("i: d= %2d o=%2d fund_d=%2d ref_o=%2d ref_d=%2d\n",
           d, octant, fundamental_d,ref_octant,ref_d) ; */
        if(R_recv_test( i_which[ref_d], &(i_plane_data[ref_d]) )
           == 0){
          printf("Null buffer not returned to DD_Sweep\n");
          error_exit(1);
        }
        for(k=0; k<i_plane_zones; k++){
          i_plane_data[ref_d][k] = eta_ref_i * i_plane_data[d][k];
        }
        i_plane_data[ref_d][i_plane_zones] = (double) ref_d;
      }

      swept[d] = 1;
      directions_left--;

    }
  }

  /* Make sure all messages have been sent */
  R_wait_send();
  return(0);
}

/*----------------------------------------------------------------------
 * CreateBufferInfoDD
 *----------------------------------------------------------------------*/

void CreateBufferInfoDD(User_Data *user_data)
{
  CreateBufferInfoDD3D(user_data);
}

/*----------------------------------------------------------------------
 * CreateBufferInfoDD3D
 *----------------------------------------------------------------------*/

void CreateBufferInfoDD3D(User_Data *user_data)
{
  Grid_Data  *grid_data  = user_data->grid_data;
  std::vector<Directions> &directions = user_data->directions;

  int *nzones          = grid_data->nzones;
  int local_imax, local_jmax, local_kmax;
  int num_directions = user_data->directions.size();
  int len[6], nm[6], length, i, d;

  local_imax = nzones[0];
  local_jmax = nzones[1];
  local_kmax = nzones[2];

  /* Info for buffers used for messages sent in the x direction */
  length = local_jmax * local_kmax + 1;
  len[0] = len[1] = length;

  /* Info for buffers used for messages sent in the y direction */
  length = local_imax * local_kmax + 1;
  len[2] = len[3] = length;

  /* Info for buffers used for messages sent in the z direction */
  length = local_imax * local_jmax + 1;
  len[4] = len[5] = length;

  for(i=0; i<6; i++){
    nm[i] = 0;
  }

  for(d=0; d<num_directions; d++){
    if(directions[d].id > 0){
      nm[0]++;
    }
    else {nm[1]++; }
    if(directions[d].jd > 0){
      nm[2]++;
    }
    else {nm[3]++; }
    if(directions[d].kd > 0){
      nm[4]++;
    }
    else {nm[5]++; }
  }

  R_buf_init( len, nm );
}
