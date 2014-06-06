/*--------------------------------------------------------------------------
 * Header file for the Directions data structures
 *--------------------------------------------------------------------------*/

#ifndef KRIPKE_DIRECTIONS_H__
#define KRIPKE_DIRECTIONS_H__

#include <vector>

class User_Data;

/**
 * Contains information needed for one quadrature set direction.
 */
struct Directions{
  double xcos;              /* Absolute value of the x-direction cosine. */
  double ycos;              /* Absolute value of the y-direction cosine. */
  double zcos;              /* Absolute value of the z-direction cosine. */
  double w;                 /* weight for the quadrature rule.*/
  int id;                   /* direction flag (= 1 if x-direction
                            cosine is positive; = -1 if not). */
  int jd;                   /* direction flag (= 1 if y-direction
                            cosine is positive; = -1 if not). */
  int kd;                   /* direction flag (= 1 if z-direction
                            cosine is positive; = -1 if not). */
  int i_src_subd;           /* Index of the neighboring spatial
                            subdomain in the negative
                            x-direction if id = 1 or positive
                            x-direction if id = -1. */
  int j_src_subd;           /* Index of the neighboring spatial
                            subdomain in the negative
                            y-direction if jd = 1 or positive
                            y-direction if jd = -1. */
  int k_src_subd;           /* Index of the neighboring spatial
                            subdomain in the negative
                            z-direction if kd = 1 or positive
                            z-direction if kd = -1. */
  int i_dst_subd;           /* Index of the neighboring spatial
                            subdomain in the positive
                            x-direction if id = 1 or negative
                            x-direction if id = -1. */
  int j_dst_subd;           /* Index of the neighboring spatial
                            subdomain in the positive
                            y-direction if jd = 1 or negative
                            y-direction if jd = -1. */
  int k_dst_subd;           /* Index of the neighboring spatial
                            subdomain in the positive
                            z-direction if kd = 1 or negative
                            z-direction if kd = -1. */
  int octant;
};


void InitDirections(User_Data *grid_data, int num_directions_per_octant);

#endif
