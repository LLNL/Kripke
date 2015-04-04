#include <Kripke/Directions.h>
#include <Kripke/Grid.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <algorithm>

/**
 * Initializes the quadrature set information for a Grid_Data object.
 * This guarantees that each <GS,DS> pair have a single originating octant.
 */
void InitDirections(Grid_Data *grid_data, int num_directions_per_octant)
{
  std::vector<Directions> &directions = grid_data->directions;

  int num_directions = 8*num_directions_per_octant;
  int omegas_per_octant=num_directions>>3;

  directions.resize(num_directions);

  for(int d=0; d<num_directions; d++){

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




