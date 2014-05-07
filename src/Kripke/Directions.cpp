#include <Kripke/Directions.h>
#include <Kripke/User_Data.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <algorithm>

static void InitOmegas(std::vector<Directions> &directions);

/* Prototypes for routines internal to this file */
struct Aux_Table {
  int index;
  int m;
  double value;
};

static bool Compare(Aux_Table const &i, Aux_Table const &j)
{
  return i.value < j.value;
}

void InitDirections(User_Data *user_data, int num_directions_per_octant)
{
  Grid_Data *grid_data = user_data->grid_data;
  std::vector<Directions> &directions = user_data->directions;
  int d, id, jd, kd, num_directions,
      direction_concurrency, full_moments;
  int in, ip, jn, jp, kn, kp;

  num_directions = 8*num_directions_per_octant;

  in = grid_data->mynbr[0][0];
  ip = grid_data->mynbr[0][1];
  jn = grid_data->mynbr[1][0];
  jp = grid_data->mynbr[1][1];
  kn = grid_data->mynbr[2][0];
  kp = grid_data->mynbr[2][1];

  directions.resize(num_directions);

  InitOmegas(directions);

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

  user_data->octant_map = octant_map;

}

void InitOmegas(std::vector<Directions> &directions)
{
  int num_directions = directions.size();
  int i, j, m, n, d, omegas_per_octant=num_directions>>3;

  std::vector<Aux_Table> aux(num_directions);

  std::vector<std::vector<double> > omegas(num_directions);
  for(m=0; m<num_directions; m++){
    omegas[m].resize(3);
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

  std::vector<int> octant_omega[8];
  int num_omegas[8];
  for(i=0; i<8; i++){
    num_omegas[i] = 0;
    octant_omega[i].resize(omegas_per_octant);
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
    std::sort(aux.begin(), aux.end(), Compare);
    for(j=0; j<omegas_per_octant; j++){
      octant_omega[i][j] = aux[j].m;
    }
  }

  std::vector<int> omega_map(num_directions);
  for(m=0, j=0; j<omegas_per_octant; j++){
    for(i=0; i<8; i++, m++){
      omega_map[m] = octant_omega[i][j];
    }
  }


  for(d=0, m=0; d<num_directions; d++, m++){

    n = omega_map[m];

    directions[d].id = (omegas[n][0] > 0.) ? 1 : -1;
    directions[d].jd = (omegas[n][1] > 0.) ? 1 : -1;
    directions[d].kd = (omegas[n][2] > 0.) ? 1 : -1;
    directions[d].xcos = fabs(omegas[n][0]);
    directions[d].ycos = fabs(omegas[n][1]);
    directions[d].zcos = fabs(omegas[n][2]);
  }

}


