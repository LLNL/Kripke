#include <Kripke/Directions.h>
#include <Kripke/User_Data.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <algorithm>

namespace {
  void InitOmegas(std::vector<Directions> &directions)
  {
    int num_directions = directions.size();
    int omegas_per_octant=num_directions>>3;

    for(int d=0, m=0; d<num_directions; d++, m++){

      int n = d/omegas_per_octant;

      double omegas[3];
      omegas[0] = n & 0x1;
      omegas[1] = (n>>1) & 0x1;
      omegas[2] = (n>>2) & 0x1;

      directions[d].id = (omegas[0] > 0.) ? 1 : -1;
      directions[d].jd = (omegas[1] > 0.) ? 1 : -1;
      directions[d].kd = (omegas[2] > 0.) ? 1 : -1;
      directions[d].xcos = fabs(omegas[0]);
      directions[d].ycos = fabs(omegas[1]);
      directions[d].zcos = fabs(omegas[2]);
    }

  }
}


/**
 * Initializes the quadrature set information for a User_Data object.
 * This guarantees that each <GS,DS> pair have a single originating octant.
 */
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

  for(d=0; d<num_directions; d++){
    id = directions[d].id;
    jd = directions[d].jd;
    kd = directions[d].kd;

    if( (id>0) && (jd>0) && (kd>0) ){
      directions[d].octant = 0;
    }
    else if( (id<0) && (jd>0) && (kd>0) ){
      directions[d].octant = 1;
    }
    else if( (id<0) && (jd<0) && (kd>0) ){
      directions[d].octant = 2;
    }
    else if( (id>0) && (jd<0) && (kd>0) ){
      directions[d].octant = 3;
    }
    else if( (id>0) && (jd>0) && (kd<0) ){
      directions[d].octant = 4;
    }
    else if( (id<0) && (jd>0) && (kd<0) ){
      directions[d].octant = 5;
    }
    else if( (id<0) && (jd<0) && (kd<0) ){
      directions[d].octant = 6;
    }
    else if( (id>0) && (jd<0) && (kd<0) ){
      directions[d].octant = 7;
    }

    directions[d].i_src_subd = (id>0) ? in : ip;
    directions[d].j_src_subd = (jd>0) ? jn : jp;
    directions[d].k_src_subd = (kd>0) ? kn : kp;
    directions[d].i_dst_subd = (id>0) ? ip : in;
    directions[d].j_dst_subd = (jd>0) ? jp : jn;
    directions[d].k_dst_subd = (kd>0) ? kp : kn;
  }
}




