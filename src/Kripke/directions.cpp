#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "transport_headers.h"
#include "comm.h"

static void InitOmegas(std::vector<Directions> &directions, int *omega_map,
                int *omega_map_inv);

/* Prototypes for routines internal to this file */
struct Aux_Table {
  int index;
  int m;
  double value;
};

static int Compare(const void *i_v, const void *j_v)
{
  struct Aux_Table *i = (struct Aux_Table *) i_v;
  struct Aux_Table *j = (struct Aux_Table *) j_v;
  /*Compare function used by qsort during omega mapping*/
  if(i->value > j->value){
    return(1);
  }
  if(i->value < j->value){
    return(-1);
  }
  return(0);
}

void InitDirections(Grid_Data *grid_data, int num_directions_per_octant)
{
  std::vector<Directions> &directions = grid_data->directions;
  int d, id, jd, kd, num_directions,
      direction_concurrency, *omega_map, *omega_map_inv,
      full_moments;
  int in, ip, jn, jp, kn, kp;

  num_directions = 8*num_directions_per_octant;

  in = grid_data->mynbr[0][0];
  ip = grid_data->mynbr[0][1];
  jn = grid_data->mynbr[1][0];
  jp = grid_data->mynbr[1][1];
  kn = grid_data->mynbr[2][0];
  kp = grid_data->mynbr[2][1];

  directions.resize(num_directions);

  NEW( omega_map, num_directions, int * );
  NEW( omega_map_inv, num_directions, int * );

  InitOmegas(directions, omega_map, omega_map_inv);

  FREE( omega_map_inv );
  FREE( omega_map );

  std::vector<int> octant_map(num_directions);

  for(d=0; d<num_directions; d++){
    id = directions[d].id;
    jd = directions[d].jd;
    kd = directions[d].kd;

    if( (id>0) && (jd>0) && (kd>0) ){
      octant_map[d] = 0;
    }
    else if( (id<0) && (jd>0) && (kd>0) ){
      octant_map[d] = 1;
    }
    else if( (id<0) && (jd<0) && (kd>0) ){
      octant_map[d] = 2;
    }
    else if( (id>0) && (jd<0) && (kd>0) ){
      octant_map[d] = 3;
    }
    else if( (id>0) && (jd>0) && (kd<0) ){
      octant_map[d] = 4;
    }
    else if( (id<0) && (jd>0) && (kd<0) ){
      octant_map[d] = 5;
    }
    else if( (id<0) && (jd<0) && (kd<0) ){
      octant_map[d] = 6;
    }
    else if( (id>0) && (jd<0) && (kd<0) ){
      octant_map[d] = 7;
    }

    directions[d].i_src_subd = (id>0) ? in : ip;
    directions[d].j_src_subd = (jd>0) ? jn : jp;
    directions[d].k_src_subd = (kd>0) ? kn : kp;
    directions[d].i_dst_subd = (id>0) ? ip : in;
    directions[d].j_dst_subd = (jd>0) ? jp : jn;
    directions[d].k_dst_subd = (kd>0) ? kp : kn;
  }

  grid_data->octant_map = octant_map;

}

void InitOmegas(std::vector<Directions> &directions, int *omega_map,
                int *omega_map_inv)
{
  int num_directions = directions.size();
  double ** omegas;
  int i, j, m, n, d, omegas_per_octant=num_directions>>3, num_omegas[8],
  *octant_omega[8];

  struct Aux_Table * aux;

  NEW( omegas, num_directions, double ** );

  NEW(aux, num_directions, struct Aux_Table * );

  for(m=0; m<num_directions; m++){
    NEW( omegas[m], 3, double * );
  }

  d = 0;
  for(i=0; i<8; i++){
    if(i == 0){
      for(m=0; m<omegas_per_octant; ++m){
        omegas[d][0] = 1.0;
        omegas[d][1] = 1.0;
        omegas[d][2] = -1.0;
        d++;
      }
    }
    else if(i == 1){
      for(m=0; m<omegas_per_octant; ++m){
        omegas[d][0] = -1.0;
        omegas[d][1] = 1.0;
        omegas[d][2] = -1.0;
        d++;
      }
    }
    else if(i == 2){
      for(m=0; m<omegas_per_octant; ++m){
        omegas[d][0] = -1.0;
        omegas[d][1] = -1.0;
        omegas[d][2] = -1.0;
        d++;
      }
    }
    else if(i == 3){
      for(m=0; m<omegas_per_octant; ++m){
        omegas[d][0] = 1.0;
        omegas[d][1] = -1.0;
        omegas[d][2] = -1.0;
        d++;
      }
    }
    else if(i == 4){
      for(m=0; m<omegas_per_octant; ++m){
        omegas[d][0] = 1.0;
        omegas[d][1] = 1.0;
        omegas[d][2] = 1.0;
        d++;
      }
    }
    else if(i == 5){
      for(m=0; m<omegas_per_octant; ++m){
        omegas[d][0] = -1.0;
        omegas[d][1] = 1.0;
        omegas[d][2] = 1.0;
        d++;
      }
    }
    else if(i == 6){
      for(m=0; m<omegas_per_octant; ++m){
        omegas[d][0] = -1.0;
        omegas[d][1] = -1.0;
        omegas[d][2] = 1.0;
        d++;
      }
    }
    else if(i == 7){
      for(m=0; m<omegas_per_octant; ++m){
        omegas[d][0] = 1.0;
        omegas[d][1] = -1.0;
        omegas[d][2] = 1.0;
        d++;
      }
    }
  }

  for(i=0; i<8; i++){
    num_omegas[i] = 0;
    NEW( octant_omega[i], omegas_per_octant, int * );
  }
/*Note that the octants are numbered in the following way (based on
  the sign of the direction cosine with the i,j, and k axis)
    i j k octant
    + + + 0        + + + 0
    - + + 1        + + - 4
    - - + 2        + - + 3
    + - + 3   ==>  + - - 7
    + + - 4        - + + 1
    - + - 5        - + - 5
    - - - 6        - - + 2
    + - - 7        - - - 6
*/
  for(m=0; m<num_directions; m++){

    if(omegas[m][0] > 0.){
      if(omegas[m][1] > 0.){
        if(omegas[m][2] > 0.){
          octant_omega[0][num_omegas[0]++] = m;
        }
        else {
          octant_omega[4][num_omegas[4]++] = m;
        }
      }
      else {
        if(omegas[m][2] > 0.){
          octant_omega[3][num_omegas[3]++] = m;
        }
        else {
          octant_omega[7][num_omegas[7]++] = m;
        }
      }
    }
    else {
      if(omegas[m][1] > 0.){
        if(omegas[m][2] > 0.){
          octant_omega[1][num_omegas[1]++] = m;
        }
        else {
          octant_omega[5][num_omegas[5]++] = m;
        }
      }
      else {
        if(omegas[m][2] > 0.){
          octant_omega[2][num_omegas[2]++] = m;
        }
        else {
          octant_omega[6][num_omegas[6]++] = m;
        }
      }
    }
  }

  /* reorder omegas in each octant using qsort */
  for(i=0; i<8; i++){
    for(j=0; j<omegas_per_octant; j++){
      aux[j].index = j;
      m = octant_omega[i][j];
      aux[j].m  = m;
      aux[j].value = 100.*fabs(omegas[m][0])
                     + 10.*fabs(omegas[m][1])
                     + fabs(omegas[m][2]);
    }
    qsort((struct Aux_Table *) aux, omegas_per_octant,
          sizeof(struct Aux_Table), Compare);
    for(j=0; j<omegas_per_octant; j++){
      octant_omega[i][j] = aux[j].m;
    }
  }
  for(m=0, j=0; j<omegas_per_octant; j++){
    for(i=0; i<8; i++, m++){
      omega_map[m] = octant_omega[i][j];
    }
  }

  for(m=0; m<num_directions; m++){
    omega_map_inv[omega_map[m]] = m;
  }

  for(d=0, m=0; d<num_directions; d++, m++){

    n = omega_map[m];

    directions[d].id = (omegas[n][0] > 0.) ? 1 : -1;
    directions[d].jd = (omegas[n][1] > 0.) ? 1 : -1;
    directions[d].kd = (omegas[n][2] > 0.) ? 1 : -1;
    directions[d].xcos = fabs(omegas[n][0]);
    directions[d].ycos = fabs(omegas[n][1]);
    directions[d].zcos = fabs(omegas[n][2]);
    /* printf("concurrency: q: %3d omega: %3d %+5.2e, %+5.2e, %+5.2e\n",
                          q,n, omegas[n][0],omegas[n][1],omegas[n][2]);  */
  }

  for(i=0; i<8; i++){
    FREE( octant_omega[i] );
  }
  for(m=0; m<num_directions; m++){
    FREE( omegas[m] );
  }
  FREE( omegas );

  FREE(aux);
}


