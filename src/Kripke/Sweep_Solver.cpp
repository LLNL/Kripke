/*--------------------------------------------------------------------------
 * Sweep-based solver routine.
 *--------------------------------------------------------------------------*/

#include <Kripke.h>
#include <Kripke/Comm.h>
#include <Kripke/User_Data.h>
#include <vector>
#include <stdio.h>



/*----------------------------------------------------------------------
 * SweepSolverSolve
 *----------------------------------------------------------------------*/

int SweepSolver (User_Data *user_data)
{
  Kernel *kernel = user_data->kernel;
  Grid_Data *grid_data = user_data->grid_data;

  BLOCK_TIMER(user_data->timing, Solve);

  // Loop over iterations
  for(int iter = 0;iter < user_data->niter;++ iter){

    /*
     * Compute the RHS
     */

    // Discrete to Moments transformation
    {
      BLOCK_TIMER(user_data->timing, LTimes);
      kernel->LTimes(grid_data);
    }


    // This is where the Scattering kernel would go!



    // Moments to Discrete transformation
    {
      BLOCK_TIMER(user_data->timing, LPlusTimes);
      kernel->LPlusTimes(grid_data);
    }

    /*
     * Sweep each Group Set
     */
    {
      BLOCK_TIMER(user_data->timing, Sweep);
      for(int group_set = 0;group_set < user_data->num_group_sets;++ group_set){
        SweepSolver_GroupSet(group_set, user_data);
      }
    }
  }
  return(0);
}


/*----------------------------------------------------------------------
 * SweepSolverSolveDD
 *----------------------------------------------------------------------*/

int SweepSolver_GroupSet (int group_set, User_Data *user_data)
{
  Grid_Data  *grid_data         = user_data->grid_data;
  Comm *comm = user_data->comm;
  std::vector<Group_Dir_Set> &dir_sets = grid_data->gd_sets[group_set];

  int num_direction_sets = user_data->num_direction_sets;

  double *msg;

  int myid;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);



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
  int emminating_directions_left;


  int in = grid_data->mynbr[0][0];
  int ip = grid_data->mynbr[0][1];
  int jn = grid_data->mynbr[1][0];
  int jp = grid_data->mynbr[1][1];
  int kn = grid_data->mynbr[2][0];
  int kp = grid_data->mynbr[2][1];

  int groups_dirs = user_data->num_groups_per_set
      * user_data->num_directions_per_set;
  int local_imax = grid_data->nzones[0];
  int local_jmax = grid_data->nzones[1];
  int local_kmax = grid_data->nzones[2];
  int i_plane_zones = local_jmax * local_kmax * groups_dirs;
  int j_plane_zones = local_imax * local_kmax * groups_dirs;
  int k_plane_zones = local_imax * local_jmax * groups_dirs;

  std::vector<double*> i_plane_data(num_direction_sets, (double*)NULL);
  std::vector<double*> j_plane_data(num_direction_sets, (double*)NULL);
  std::vector<double*> k_plane_data(num_direction_sets, (double*)NULL);

  std::vector<int> i_which(num_direction_sets);
  std::vector<int> j_which(num_direction_sets);
  std::vector<int> k_which(num_direction_sets);
  for(int ds=0; ds<num_direction_sets; ds++){
    Directions *directions = dir_sets[ds].directions;
    i_which[ds] = (directions[0].id>0) ? 0 : 1;
    j_which[ds] = (directions[0].jd>0) ? 2 : 3;
    k_which[ds] = (directions[0].kd>0) ? 4 : 5;
  }

  std::vector<double> psi_lf_data((local_imax+1)*local_jmax*local_kmax);
  std::vector<double> psi_fr_data( local_imax*(local_jmax+1)*local_kmax);
  std::vector<double> psi_bo_data(local_imax*local_jmax*(local_kmax+1));

  /* Hang out receive requests for each of the 6 neighbors */
  if(in != -1){
    comm->R_recv_dir( 0, in );
  }
  if(ip != -1){
    comm->R_recv_dir( 1, ip );
  }
  if(jn != -1){
    comm->R_recv_dir( 2, jn );
  }
  if(jp != -1){
    comm->R_recv_dir( 3, jp );
  }
  if(kn != -1){
    comm->R_recv_dir( 4, kn );
  }
  if(kp != -1){
    comm->R_recv_dir( 5, kp );
  }

  /* Allocate and initialize (set to zero for now) message
     buffers for subdomain faces on the problem boundary */
  for(int ds=0; ds<num_direction_sets; ds++){

    Directions *directions = dir_sets[ds].directions;
    int i_src_subd = directions[0].i_src_subd;
    int j_src_subd = directions[0].j_src_subd;
    int k_src_subd = directions[0].k_src_subd;

    if(k_src_subd == -1){
      if(comm->R_recv_test( k_which[ds], &(k_plane_data[ds]) ) == 0){
        printf("Null buffer not returned to DD_Sweep\n");
        error_exit(1);
      }
      for(int k=0; k<k_plane_zones; k++){
        k_plane_data[ds][k] = 1.0;
      }
      k_plane_data[ds][k_plane_zones] = (double) ds;
    }
    else {
      k_plane_data[ds] = NULL;
    }

    if(j_src_subd == -1){
      if(comm->R_recv_test( j_which[ds], &(j_plane_data[ds]) ) == 0){
        printf("Null buffer not returned to DD_Sweep\n");
        error_exit(1);
      }
      for(int k=0; k<j_plane_zones; k++){
        j_plane_data[ds][k] = 1.0;
      }
      j_plane_data[ds][j_plane_zones] = (double) ds;
    }
    else {
      j_plane_data[ds] = NULL;
    }

    if(i_src_subd == -1){
      if(comm->R_recv_test( i_which[ds], &(i_plane_data[ds]) ) == 0){
        printf("Null buffer not returned to DD_Sweep\n");
        error_exit(1);
      }
      for(int i=0; i<i_plane_zones; i++){
        i_plane_data[ds][i] = 1.0;
      }
      i_plane_data[ds][i_plane_zones] = (double) ds;
    }
    else {
      i_plane_data[ds] = NULL;
    }
  }

  int directions_left = num_direction_sets;
  std::vector<int> swept(num_direction_sets, 0.0);

  while(directions_left){

    /* Check for a message from the 6 neighboring subdomains. */
    if(in != -1 && comm->R_recv_test( 0, &msg ) != 0){
      i_plane_data[(int) msg[i_plane_zones]] = msg;
    }

    if(ip != -1 && comm->R_recv_test( 1, &msg ) != 0){
      i_plane_data[(int) msg[i_plane_zones]] = msg;
    }

    if(jn != -1 && comm->R_recv_test( 2, &msg ) != 0){
      j_plane_data[(int) msg[j_plane_zones]] = msg;
    }

    if(jp != -1 && comm->R_recv_test( 3, &msg ) != 0){
      j_plane_data[(int) msg[j_plane_zones]] = msg;
    }

    if(kn != -1 && comm->R_recv_test( 4, &msg ) != 0){
      k_plane_data[(int) msg[k_plane_zones]] = msg;
    }

    if(kp != -1 && comm->R_recv_test( 5, &msg ) != 0){
      k_plane_data[(int) msg[k_plane_zones]] = msg;
    }

    for(int ds=0; ds<num_direction_sets; ds++){
      if(k_plane_data[ds] == NULL ||
         j_plane_data[ds] == NULL ||
         i_plane_data[ds] == NULL ||
         swept[ds]){
        continue;
      }

      /* Use standard Diamond-Difference sweep */
      {
        BLOCK_TIMER(user_data->timing, Sweep_Kernel);
        user_data->kernel->sweep(grid_data, &dir_sets[ds], i_plane_data[ds], j_plane_data[ds], k_plane_data[ds]);
      }

      Directions *directions = dir_sets[ds].directions;
      int i_dst_subd = directions[0].i_dst_subd;
      int j_dst_subd = directions[0].j_dst_subd;
      int k_dst_subd = directions[0].k_dst_subd;

      comm->R_send( i_plane_data[ds], i_dst_subd, i_plane_zones+1 );
      comm->R_send( j_plane_data[ds], j_dst_subd, j_plane_zones+1 );
      comm->R_send( k_plane_data[ds], k_dst_subd, k_plane_zones+1 );

      swept[ds] = 1;
      directions_left--;

    }
  }

  /* Make sure all messages have been sent */
  comm->R_wait_send();
  return(0);
}

/*----------------------------------------------------------------------
 * CreateBufferInfo
 *----------------------------------------------------------------------*/

void CreateBufferInfo(User_Data *user_data)
{
  Grid_Data  *grid_data  = user_data->grid_data;

  int *nzones          = grid_data->nzones;
  int local_imax, local_jmax, local_kmax;
  int num_direction_sets = grid_data->gd_sets[0].size();
  int len[6], nm[6], length;

  // get group and direction dimensionality
  int dirs_groups = user_data->num_directions_per_set
                  * user_data->num_groups_per_set;

  local_imax = nzones[0];
  local_jmax = nzones[1];
  local_kmax = nzones[2];

  /* Info for buffers used for messages sent in the x direction */
  length = local_jmax * local_kmax + 1;
  len[0] = len[1] = length * dirs_groups;

  /* Info for buffers used for messages sent in the y direction */
  length = local_imax * local_kmax + 1;
  len[2] = len[3] = length * dirs_groups;

  /* Info for buffers used for messages sent in the z direction */
  length = local_imax * local_jmax + 1;
  len[4] = len[5] = length * dirs_groups;

  for(int i=0; i<6; i++){
    nm[i] = 0;
  }

  for(int ds=0; ds<num_direction_sets; ds++){
    Directions *directions = grid_data->gd_sets[0][ds].directions;
    if(directions[0].id > 0){
      nm[0]++;
    }
    else {nm[1]++; }
    if(directions[0].jd > 0){
      nm[2]++;
    }
    else {nm[3]++; }
    if(directions[0].kd > 0){
      nm[4]++;
    }
    else {nm[5]++; }
  }

  user_data->comm = new Comm( len, nm );
}
