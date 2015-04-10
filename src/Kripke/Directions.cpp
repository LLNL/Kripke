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
  directions.resize(num_directions);

  int d = 0;
  for(int octant = 0;octant < 8;++ octant){
    double omegas[3];
    omegas[0] = octant & 0x1;
    omegas[1] = (octant>>1) & 0x1;
    omegas[2] = (octant>>2) & 0x1;

    for(int sd=0; sd<num_directions_per_octant; sd++, d++){
      // Store which logical direction of travel we have
      directions[d].id = (omegas[0] > 0.) ? 1 : -1;
      directions[d].jd = (omegas[1] > 0.) ? 1 : -1;
      directions[d].kd = (omegas[2] > 0.) ? 1 : -1;


      // Store quadrature point's weight
      directions[d].w = 1.0 / (double)num_directions;

      // Get the direction of this quadrature point

      double theta = M_PI/4;
      double omega = M_PI/4;
      /*
      double theta, omega;
      switch(sd){
      case 0:
          theta = acos(0.3500212);
          omega = M_PI/4;
          break;
      case 1:
          omega = (M_PI/2)/3;
      case 2:
          omega = 2*(M_PI/2)/3;
          theta = acos(0.8688903);
      }*/


      // Compute x,y,z cosine values
      double mu  = cos(theta);
      double eta = sqrt(1-mu*mu) * cos(omega);
      double xi  = sqrt(1-mu*mu) * sin(omega);
      directions[d].xcos = mu;
      directions[d].ycos = eta;
      directions[d].zcos = xi;
    }
  }
}




