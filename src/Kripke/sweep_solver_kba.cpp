/*--------------------------------------------------------------------------
 * Sweep-based solver routine.
 *--------------------------------------------------------------------------*/

#include "transport_headers.h"

/*----------------------------------------------------------------------
 * SweepSolverSolveDDKBA
 *----------------------------------------------------------------------*/

int SweepSolverSolveDDKBA(User_Data *user_data, double **rhs, double **ans,
			  double **tempv)
{
  Grid_Data          *grid_data         = user_data->grid_data;
  Grid_Data          *grid_data_kba;
  Directions         *directions        = grid_data->directions;
  Boltzmann_Solver  *boltzmann_solver = user_data->boltzmann_solver;
  Sweep_Solver_Data *sweep_solver_data =
    boltzmann_solver->sweep_solver_data;
  Sigma_Tot          *sigma_tot         = user_data->sigma_tot;
  Data_Vector        *tmp_sigma_tot     = user_data->tmp_sigma_tot;

  double    *volume            = grid_data->volume;
  double   **psi_i_plane       = user_data->psi_i_plane;
  double   **psi_j_plane       = user_data->psi_j_plane;
  double   **psi_k_plane       = user_data->psi_k_plane;

  int *nzones = grid_data->nzones;
  int nlevels_kba = user_data->nlevels_kba;
  int num_zones = nzones[0]*nzones[1]*nzones[2];
  int nx=nzones[0], ny=nzones[1], nz=nzones[2];
  int num_directions = grid_data->num_directions;

  double *msg, **i_plane_data, **j_plane_data, **k_plane_data;
  int i, k, sweep_group, *i_plane_zones, *j_plane_zones, k_plane_zones,
      i_src_subd, j_src_subd, k_src_subd, i_dst_subd,
      j_dst_subd, k_dst_subd, directions_and_levels_left, in, ip,
      jn, jp, kn, kp, d, *i_which, *j_which, *k_which;
  int local_imax, local_jmax, *local_kmax;
  int *offset, *dz_offset, *i_plane_zones_offset, *j_plane_zones_offset;
  int ref_dl, dl, kl, lev;

  int *swept; // Need local array, not sweep_solver_data->swept;
  int bc_ref_in, bc_ref_ip, bc_ref_jn, bc_ref_jp, bc_ref_kn, bc_ref_kp;
  double eta_ref_in, eta_ref_ip, eta_ref_jn, eta_ref_jp;
  double eta_ref_kn, eta_ref_kp;
  double eta_ref_i, eta_ref_j, eta_ref_k;

  double *psi_lf_data, *psi_fr_data, *psi_bo_data;

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
  int emminating_directions_and_levels_left;

  NEW(swept, num_directions*nlevels_kba, int *);
  NEW(grid_data_kba, 1, Grid_Data *); /* Also need nzones here */
  NEW(grid_data_kba->nzones, 3, int *);

  bc_ref_in = user_data->bc_data->bc_types[0];
  bc_ref_ip = user_data->bc_data->bc_types[1];
  bc_ref_jn = user_data->bc_data->bc_types[2];
  bc_ref_jp = user_data->bc_data->bc_types[3];
  bc_ref_kn = user_data->bc_data->bc_types[4];
  bc_ref_kp = user_data->bc_data->bc_types[5];

  in = grid_data->mynbr[0][0];
  ip = grid_data->mynbr[0][1];
  jn = grid_data->mynbr[1][0];
  jp = grid_data->mynbr[1][1];
  kn = grid_data->mynbr[2][0];
  kp = grid_data->mynbr[2][1];

  if(in == -1){
    eta_ref_in = user_data->bc_data->bc_values[0];
  }
  else {
    eta_ref_in = 0.0;
  }
  if(ip == -1){
    eta_ref_ip = user_data->bc_data->bc_values[1];
  }
  else {
    eta_ref_ip = 0.0;
  }
  if(jn == -1){
    eta_ref_jn = user_data->bc_data->bc_values[2];
  }
  else {
    eta_ref_jn = 0.0;
  }
  if(jp == -1){
    eta_ref_jp = user_data->bc_data->bc_values[3];
  }
  else {
    eta_ref_jp = 0.0;
  }
  if(kn == -1){
    eta_ref_kn = user_data->bc_data->bc_values[4];
  }
  else {
    eta_ref_kn = 0.0;
  }
  if(kp == -1){
    eta_ref_kp = user_data->bc_data->bc_values[5];
  }
  else {
    eta_ref_kp = 0.0;
  }

  local_imax = nzones[0];
  local_jmax = nzones[1];

  /* Calculate local_kmax, i_plane_zones, j_plane_zones, and k_plane_zones.
     Note that local_kmax[0] contains the maximum value of all the
     local_kmax's. Also, i_planes_zones[0] and j_plane_zones[0] contain
     the maximum of all the i_plane_zones and j_plane_zones, respectively. */
  NEW(local_kmax, nlevels_kba, int *);
  NEW(i_plane_zones, nlevels_kba, int*);
  NEW(j_plane_zones, nlevels_kba, int*);
  int base = nzones[2]/nlevels_kba;
  int rem = nzones[2] % nlevels_kba;
  for(i = 0; i < nlevels_kba; i++){
    local_kmax[i] = base;
    if(i < rem){
      local_kmax[i] += 1;
    }
    i_plane_zones[i] = local_jmax * local_kmax[i];
    j_plane_zones[i] = local_imax * local_kmax[i];
  }
  k_plane_zones = local_imax * local_jmax;

  /* Calculate zonal offset for each level */
  NEW(offset, nlevels_kba, int *);
  offset[0] = 0;
  for(lev=1; lev<nlevels_kba; lev++){
    offset[lev] = offset[lev-1] + local_imax * local_jmax * local_kmax[lev-1];
  }

  /* Calculate dz_offset for each level */
  NEW(dz_offset, nlevels_kba, int *);
  dz_offset[0] = 0;
  for(lev=1; lev<nlevels_kba; lev++){
    dz_offset[lev] = dz_offset[lev-1] + local_kmax[lev-1];
  }

  /* Calculate i_plane_zones_offset and j_plane_zones_offset for
     each level */
  NEW(i_plane_zones_offset, nlevels_kba, int *);
  NEW(j_plane_zones_offset, nlevels_kba, int *);
  i_plane_zones_offset[0] = 0;
  j_plane_zones_offset[0] = 0;
  for(lev=1; lev<nlevels_kba; lev++){
    i_plane_zones_offset[lev] = i_plane_zones_offset[lev-1] +
                                i_plane_zones[lev-1];
    j_plane_zones_offset[lev] = j_plane_zones_offset[lev-1] +
                                j_plane_zones[lev-1];
  }

  NEW( i_plane_data, num_directions*nlevels_kba, double ** );
  NEW( j_plane_data, num_directions*nlevels_kba, double ** );
  NEW( k_plane_data, num_directions*nlevels_kba, double ** );

  NEW( i_which, num_directions*nlevels_kba, int * );
  NEW( j_which, num_directions*nlevels_kba, int * );
  NEW( k_which, num_directions*nlevels_kba, int * );

  for(lev=0; lev<nlevels_kba; lev++){
    for(d=0; d<num_directions; d++){
      i_which[d+lev*num_directions] = (directions[d].id>0) ? 0 : 1;
      j_which[d+lev*num_directions] = (directions[d].jd>0) ? 2 : 3;
      k_which[d+lev*num_directions] = (directions[d].kd>0) ? 4 : 5;
    }
  }

  /* Need to use maximum local_kmax here */
  NEW( psi_lf_data, (local_imax+1)*local_jmax*local_kmax[0], double * );
  NEW( psi_fr_data, local_imax*(local_jmax+1)*local_kmax[0], double * );
  NEW( psi_bo_data, local_imax*local_jmax*(local_kmax[0]+1), double * );

  /* Evaluate the total cross section for this group */
  EvalSigmaTot(sigma_tot, tmp_sigma_tot);

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
  for(lev=0; lev<nlevels_kba; lev++){
    for(d=0; d<num_directions; d++){
      dl = d + lev*num_directions; /* Offset into direction and level */

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

      if(k_src_subd == -1 && bc_ref_k == 0 &&
         (((lev == 0) && (directions[d].kd > 0)) ||
          ((lev == nlevels_kba-1) && (directions[d].kd < 0)))){
        if(R_recv_test( k_which[dl], &(k_plane_data[dl]) ) == 0){
          printf("Null buffer not returned to DD_Sweep\n");
          error_exit(1);
        }
        for(k=0; k<k_plane_zones; k++){
          k_plane_data[dl][k] = eta_ref_k;
        }
        k_plane_data[dl][k_plane_zones] = (double) dl;
      }
      else {
        k_plane_data[dl] = NULL;
      }

      if(j_src_subd == -1 && bc_ref_j == 0){
        if(R_recv_test( j_which[dl], &(j_plane_data[dl]) ) == 0){
          printf("Null buffer not returned to DD_Sweep\n");
          error_exit(1);
        }
        for(k=0; k<j_plane_zones[lev]; k++){
          j_plane_data[dl][k] = eta_ref_j;
        }
        j_plane_data[dl][j_plane_zones[0]] = (double) dl;
      }
      else {
        j_plane_data[dl] = NULL;
      }

      if(i_src_subd == -1 && bc_ref_i == 0){
        if(R_recv_test( i_which[dl], &(i_plane_data[dl]) ) == 0){
          printf("Null buffer not returned to DD_Sweep\n");
          error_exit(1);
        }
        for(i=0; i<i_plane_zones[lev]; i++){
          i_plane_data[dl][i] = eta_ref_i;
        }
        i_plane_data[dl][i_plane_zones[0]] = (double) dl;
      }
      else {
        i_plane_data[dl] = NULL;
      }
    }
  }

  directions_and_levels_left = num_directions*nlevels_kba;
  for(dl=0; dl<num_directions*nlevels_kba; dl++){
    swept[dl] = 0;
  }
  while(directions_and_levels_left){

    /* Check for a message from the 6 neighboring subdomains. */
    if(in != -1 && R_recv_test( 0, &msg ) != 0){
      i_plane_data[(int) msg[i_plane_zones[0]]] = msg;
    }

    if(ip != -1 && R_recv_test( 1, &msg ) != 0){
      i_plane_data[(int) msg[i_plane_zones[0]]] = msg;
    }

    if(jn != -1 && R_recv_test( 2, &msg ) != 0){
      j_plane_data[(int) msg[j_plane_zones[0]]] = msg;
    }

    if(jp != -1 && R_recv_test( 3, &msg ) != 0){
      j_plane_data[(int) msg[j_plane_zones[0]]] = msg;
    }

    if(kn != -1 && R_recv_test( 4, &msg ) != 0){
      k_plane_data[(int) msg[k_plane_zones]] = msg;
    }

    if(kp != -1 && R_recv_test( 5, &msg ) != 0){
      k_plane_data[(int) msg[k_plane_zones]] = msg;
    }

    for(lev=0; lev<nlevels_kba; lev++){
      for(d=0; d<num_directions; d++){
        dl = d + lev*num_directions; /* Offset into direction and level
                                             */
        if(k_plane_data[dl] == NULL ||
           j_plane_data[dl] == NULL ||
           i_plane_data[dl] == NULL ||
           swept[dl]){
          continue;
        }

        /* Need to load a temporary grid_data structure */
        grid_data_kba->nzones[0] = local_imax;
        grid_data_kba->nzones[1] = local_jmax;
        grid_data_kba->nzones[2] = local_kmax[lev];

        grid_data_kba->directions = grid_data->directions;

        grid_data_kba->deltas[0] = grid_data->deltas[0];
        grid_data_kba->deltas[1] = grid_data->deltas[1];
        grid_data_kba->deltas[2] = grid_data->deltas[2] + dz_offset[lev];

	/* Use standard Diamond-Difference sweep */
	SweepDD(d, grid_data_kba, volume+offset[lev],
		tmp_sigma_tot->data+offset[lev],
		rhs[d]+offset[lev], ans[d]+offset[lev],
		i_plane_data[dl], j_plane_data[dl],
		k_plane_data[dl],
		psi_lf_data, psi_fr_data, psi_bo_data);

        /* k_plane_data is copied for next level sweep */
        if((directions[d].kd > 0) && (lev < nlevels_kba-1)){
          int dlp1 = d + (lev+1)*num_directions;
          if(R_recv_test( k_which[dlp1], &(k_plane_data[dlp1]) ) == 0){
            printf("Null buffer not returned to DD_Sweep\n");
            error_exit(1);
          }
          /* printf("dlp1=%d\n",dlp1); */
          for(k=0; k<k_plane_zones; k++){
            k_plane_data[dlp1][k] = k_plane_data[dl][k];
          }
          k_plane_data[dlp1][k_plane_zones] = (double) dlp1;
        }
        else if((directions[d].kd < 0) && (lev > 0)){
          int dlm1 = d + (lev-1)*num_directions;
          if(R_recv_test( k_which[dlm1], &(k_plane_data[dlm1]) ) == 0){
            printf("Null buffer not returned to DD_Sweep\n");
            error_exit(1);
          }
          /* printf("dlm1=%d\n",dlm1); */
          for(k=0; k<k_plane_zones; k++){
            k_plane_data[dlm1][k] = k_plane_data[dl][k];
          }
          k_plane_data[dlm1][k_plane_zones] = (double) dlm1;
        }

        i_dst_subd = directions[d].i_dst_subd;
        j_dst_subd = directions[d].j_dst_subd;
        k_dst_subd = directions[d].k_dst_subd;

        R_send( i_plane_data[dl], i_dst_subd, i_plane_zones[0]+1 );
        R_send( j_plane_data[dl], j_dst_subd, j_plane_zones[0]+1 );
        R_send( k_plane_data[dl], k_dst_subd, k_plane_zones+1 );

        /*Check if any of the 3 planes are reflective problem boundaries.
          If so, generate the src for future sweeps */

        /* get reflective b.c. information for dst faces */
        bc_ref_i = (directions[d].id>0) ? bc_ref_ip : bc_ref_in;
        bc_ref_j = (directions[d].jd>0) ? bc_ref_jp : bc_ref_jn;
        bc_ref_k = (directions[d].kd>0) ? bc_ref_kp : bc_ref_kn;
        eta_ref_i = (directions[d].id>0) ? eta_ref_ip : eta_ref_in;
        eta_ref_j = (directions[d].jd>0) ? eta_ref_jp : eta_ref_jn;
        eta_ref_k = (directions[d].kd>0) ? eta_ref_kp : eta_ref_kn;

        if(k_dst_subd == -1 && bc_ref_k == 1 &&
           ((lev == 0 && bc_ref_k == bc_ref_kn) ||
            (lev == nlevels_kba-1 && bc_ref_k == bc_ref_kp)) ){
          octant = directions->octant_map[d];
          ref_octant = r_rules[octant][2];
	  fundamental_d = (d - octant)/8;
	  ref_d = 8 * fundamental_d + ref_octant;
          /* printf("k: d= %2d o=%2d fund_d=%2d ref_o=%2d ref_d=%2d\n",
             d, octant, fundamental_d,ref_octant,ref_d) ; */
          ref_dl = ref_d + lev*num_directions;
          if(R_recv_test( k_which[ref_dl], &(k_plane_data[ref_dl]) )
             == 0){
            printf("Null buffer not returned to DD_Sweep\n");
            error_exit(1);
          }
          for(k=0; k<k_plane_zones; k++){
            k_plane_data[ref_dl][k] = eta_ref_k * k_plane_data[dl][k];
          }
          k_plane_data[ref_dl][k_plane_zones] = (double) ref_dl;
        }
        if(j_dst_subd == -1 && bc_ref_j == 1){
          octant = directions->octant_map[d];
          ref_octant = r_rules[octant][1];
	  fundamental_d = (d - octant)/8;
	  ref_d = 8 * fundamental_d + ref_octant;
          /* printf("j: d= %2d o=%2d fund_d=%2d ref_o=%2d ref_d=%2d\n",
             d, octant, fundamental_d,ref_octant,ref_d) ; */
          ref_dl = ref_d + lev*num_directions;
          if(R_recv_test( j_which[ref_dl], &(j_plane_data[ref_dl]) )
             == 0){
            printf("Null buffer not returned to DD_Sweep\n");
            error_exit(1);
          }
          for(k=0; k<j_plane_zones[lev]; k++){
            j_plane_data[ref_dl][k] = eta_ref_j * j_plane_data[dl][k];
          }
          j_plane_data[ref_dl][j_plane_zones[0]] = (double) ref_dl;
        }
        if(i_dst_subd == -1 && bc_ref_i == 1){
          octant = directions->octant_map[d];
          ref_octant = r_rules[octant][0];
	  fundamental_d = (d - octant)/8;
	  ref_d = 8 * fundamental_d + ref_octant;
          /* printf("i: d= %2d o=%2d fund_d=%2d ref_o=%2d ref_d=%2d\n",
             d, octant, fundamental_d,ref_octant,ref_d) ; */
          ref_dl = ref_d + lev*num_directions;
          if(R_recv_test( i_which[ref_dl], &(i_plane_data[ref_dl]) )
             == 0){
            printf("Null buffer not returned to DD_Sweep\n");
            error_exit(1);
          }
          for(k=0; k<i_plane_zones[lev]; k++){
            i_plane_data[ref_dl][k] = eta_ref_i * i_plane_data[dl][k];
          }
          i_plane_data[ref_dl][i_plane_zones[0]] = (double) ref_dl;
        }

        swept[dl] = 1;
        directions_and_levels_left--;

        /* Copy outgoing face data into psi_i_plane, psi_j_plane and
           psi_k_plane. These are needed for leakage calculation. */
        for(k=0; k<i_plane_zones[lev]; k++){
          kl = k + i_plane_zones_offset[lev];
          psi_i_plane[d][kl] = i_plane_data[dl][k];
        }
        for(k=0; k<j_plane_zones[lev]; k++){
          kl = k + j_plane_zones_offset[lev];
          psi_j_plane[d][kl] = j_plane_data[dl][k];
        }
        if((lev == 0) || (lev == nlevels_kba-1)){
          for(k=0; k<k_plane_zones; k++){
            psi_k_plane[d][k] = k_plane_data[dl][k];
          }
        }

      } /* End of loop over directions */
    } /* End of loop over kba levels */
  } /* End of while loop */

  /* Make sure all messages have been sent */
  R_wait_send();

  /* SynchronizeR(); */

  /* Free temporary space */
  FREE( i_which );
  FREE( j_which );
  FREE( k_which );
  FREE( i_plane_data );
  FREE( j_plane_data );
  FREE( k_plane_data );

  FREE( psi_lf_data );
  FREE( psi_fr_data );
  FREE( psi_bo_data );

  FREE( local_kmax );
  FREE( offset );
  FREE( i_plane_zones );
  FREE( j_plane_zones );
  FREE( i_plane_zones_offset );
  FREE( j_plane_zones_offset );
  FREE( dz_offset );

  FREE(swept);
  FREE(grid_data_kba->nzones);
  FREE(grid_data_kba);

  return(0);
}

/*----------------------------------------------------------------------
 * CreateBufferInfoDDKBA
 *----------------------------------------------------------------------*/

void CreateBufferInfoDDKBA(User_Data *user_data)
{
  Grid_Data *grid_data  = user_data->grid_data;
  Directions *directions = grid_data->directions;

  int *nzones = grid_data->nzones;
  int local_imax, local_jmax;
  int num_directions = grid_data->num_directions;
  int nlevels_kba = user_data->nlevels_kba;
  int len[6], nm[6], length, i, d, lev;
  int *local_kmax;
  int base, rem;

  NEW(local_kmax, nlevels_kba, int *);

  local_imax = nzones[0];
  local_jmax = nzones[1];

  base = nzones[2]/nlevels_kba;
  rem = nzones[2] % nlevels_kba;

  for(i = 0; i < nlevels_kba; i++){
    local_kmax[i] = base;
    if(i < rem){
      local_kmax[i] += 1;
    }
  }

  /* Info for buffers used for messages sent in the x direction */
  length = local_jmax * local_kmax[0] + 1;
  len[0] = len[1] = length;

  /* Info for buffers used for messages sent in the y direction */
  length = local_imax * local_kmax[0] + 1;
  len[2] = len[3] = length;

  /* Info for buffers used for messages sent in the z direction */
  length = local_imax * local_jmax + 1;
  len[4] = len[5] = length;

  for(i=0; i<6; i++){
    nm[i] = 0;
  }

  for(lev=0; lev<nlevels_kba; lev++){
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
  }

  R_buf_init( len, nm );

  FREE( local_kmax );
}
