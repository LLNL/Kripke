/*
 * NOTICE
 *
 * This work was produced at the Lawrence Livermore National Laboratory (LLNL)
 * under contract no. DE-AC-52-07NA27344 (Contract 44) between the U.S.
 * Department of Energy (DOE) and Lawrence Livermore National Security, LLC
 * (LLNS) for the operation of LLNL. The rights of the Federal Government are
 * reserved under Contract 44.
 *
 * DISCLAIMER
 *
 * This work was prepared as an account of work sponsored by an agency of the
 * United States Government. Neither the United States Government nor Lawrence
 * Livermore National Security, LLC nor any of their employees, makes any
 * warranty, express or implied, or assumes any liability or responsibility
 * for the accuracy, completeness, or usefulness of any information, apparatus,
 * product, or process disclosed, or represents that its use would not infringe
 * privately-owned rights. Reference herein to any specific commercial products,
 * process, or service by trade name, trademark, manufacturer or otherwise does
 * not necessarily constitute or imply its endorsement, recommendation, or
 * favoring by the United States Government or Lawrence Livermore National
 * Security, LLC. The views and opinions of authors expressed herein do not
 * necessarily state or reflect those of the United States Government or
 * Lawrence Livermore National Security, LLC, and shall not be used for
 * advertising or product endorsement purposes.
 *
 * NOTIFICATION OF COMMERCIAL USE
 *
 * Commercialization of this product is prohibited without notifying the
 * Department of Energy (DOE) or Lawrence Livermore National Security.
 */

#ifndef __CLIPLIANGBARSKY3D_HH__
#define __CLIPLIANGBARSKY3D_HH__

#include <cmath>
#include <limits>
#include <array>

namespace KKC_CLIP_3D {

  // Implmentation of Liang Barsky algorithm for clipping against coordinate aligned Cartesian boxes
  // This code essentially follows that in :
  //   Foley, Van Dam, Feiner, and Hughes, "Computer Graphics", 2nd Edition in C, 1996.
  /* We make a slight modification in an attempt to make the clipping more robust to degenerate intersections and lines parallel to edges.
     Basically, we "fuzz" up the box near machine epsilon (relative to the scale of the box coordinates). The box can either be a bit bigger than it
     really is or a bit smaller depending on the sign of a template argument.  Note this also affects whether the box boundary closes the box volume or not.
     If we were really worried about this, we would use either robust predicates (adaptive precision math) 
     or an implementation of Simulation of Simplicity. The implementation here is basically using tolerances.
  */

  // this is just a helper to clean things up a bit in the code, there is an interface that avoids calling codes needing to know about this
  template<typename T>
  using Vector3D = std::array<T,3>;
  
  template<typename T>
  inline
  auto
  clip_or_reject(const T &norm_dot_line, // input: the cutting plane's outward normal dotted with the line
		 const T &norm_dot_dir,  // input: the same outward normal dotted with a line going from the start point to a point on the cutting plane
		 const T &line_scale, // input: the scaling used in the tolerance of the normal dotted with the line
		 const T &dir_scale,  // input: the scaling used in the tolerance of the normal dotted with the direction vector
		 T &t_min, T &t_max) // input and output: the parameter space endpoints of the clipped line, if it is clipped
  {
    // return output: true if we can accept this line or false if we might be able to reject it
    bool accept = true;
    T t;

    T denom = -norm_dot_line; // I like making the input to this function clear, so do this to make the logic follow Foley et al.
    if ( denom > std::abs(line_scale) ) 
      { // entering itersection
	t = norm_dot_dir/denom;
	if ( (t-t_max)>std::numeric_limits<T>::epsilon() ) // note that t\in[0,1] is usually an O(1) number, so this comparison is appropriate, however we should really use a robust predicate 
	  accept=false;
	else if ( (t-t_min)>std::numeric_limits<T>::epsilon() )
	  t_min = t; // found a new t_min on the line segment
      }
    else if ( denom < -std::abs(line_scale) ) 
      { // leaving intersection
	t = norm_dot_dir/denom;
	if ( (t_min-t)>std::numeric_limits<T>::epsilon() ) // note that t\in[0,1] is usually an O(1) number, so this comparison is appropriate, however we should really use a robust predicate 
	  accept=false;
	else if ( (t_max-t)>std::numeric_limits<T>::epsilon() )
	  t_max = t; // found a new t_max on the line segment
      }
    else if ( norm_dot_dir > dir_scale ) 
      { 
	accept = false;
      }
    return accept;
  }
  
  template<typename T, int REL_EPS=5> // increase REL_EPS to make the inersection calculation more "sloppy", note making REL_EPS negative effectively shrinks the box
  inline
  auto
  LiangBarskyClip3DBox(const Vector3D<T> &box_lower_bounds, const Vector3D<T> &box_upper_bounds, // input: bounding points for the clipping box
		       const Vector3D<T> &p0, const Vector3D<T> &p1, // input: the end points of the line segment to clip
		       T& t_min, T&t_max) // output: the line paramter coordinates of the clipping endpoints 
  {
    // return output: true if the line has a segment inside the box, false if not
    bool segment_inside_box = false;
    
    t_min = 0.0;
    t_max = 1.0;

    Vector3D<T> dx, dx_box;
    T unperturbed_l1_norm=0, coord_norm=0, box_norm=0; // norms used for various tolerances
    for ( size_t i=0; i<3; i++ )
      {
	unperturbed_l1_norm += std::abs(p1[i] - p0[i]);

	// strictly speaking we don't need the following stuff yet, but lets compute it anyways
	coord_norm += (std::abs(p1[i])+std::abs(p0[i]))/6.0; // this will give us an estimate for the average scale of the coordinates
	box_norm += (std::abs(box_lower_bounds[i]) + std::abs(box_upper_bounds[i]) + p0[i])/9.0;
	dx[i] = p1[i] - p0[i];
	dx_box[i] = box_upper_bounds[i] - box_lower_bounds[i];
      }
    
    // check to see if the line segment is degenerate
    if (unperturbed_l1_norm < std::numeric_limits<T>::epsilon()*coord_norm || unperturbed_l1_norm<std::numeric_limits<T>::min() )
      { // this line segment is degenerate, really a point
	//   check to see if this thing is in or out of the box, use only one of the unperturbed points
	bool is_inside_box = true;
	for ( size_t i=0; i<3 && is_inside_box; i++ )
	  is_inside_box = (p0[i]>=box_lower_bounds[i]) && (p0[i]<=box_upper_bounds[i]); //yeah, sigh, should use scaled epsilon here too...

	// now get the heck out of here
	if (is_inside_box)
	  return true;
	else
	  return false;
      }

    const T line_scale = coord_norm * REL_EPS * std::numeric_limits<T>::epsilon();
    const T box_scale  = box_norm * REL_EPS * std::numeric_limits<T>::epsilon();
    // continue with the line segment clipping since we now know we have an actual line segment
    // note that for the left, bottom and back faces we use the box lower bounds for the direction calculation
    //           for the right, top and front faces we use the box upper bounds for the direction calculation
    if (clip_or_reject( -dx[0], -(p0[0]-box_lower_bounds[0]), line_scale, box_scale, t_min, t_max) ) // left face
      if (clip_or_reject(  dx[0], (p0[0]-box_upper_bounds[0]), line_scale, box_scale, t_min, t_max) ) // right face
	if (clip_or_reject( -dx[1], -(p0[1]-box_lower_bounds[1]), line_scale, box_scale, t_min, t_max) ) // bottom face
	  if (clip_or_reject( dx[1], (p0[1]-box_upper_bounds[1]), line_scale, box_scale, t_min, t_max) ) // top face
	    if (clip_or_reject( -dx[2], -(p0[2]-box_lower_bounds[2]), line_scale, box_scale, t_min, t_max) ) // back face
	      if (clip_or_reject(  dx[2], (p0[2]-box_upper_bounds[2]), line_scale, box_scale, t_min, t_max) ) // front face
		segment_inside_box = std::abs(t_max-t_min)>std::numeric_limits<T>::epsilon(); // Foley et al don't check this but I do...

    return segment_inside_box;
  }

  // here is an interface that does not use "internal" data structures like Vector3D
  
  template<typename T, int REL_EPS=2> // increase REL_EPS to make the inersection calculation more "sloppy", note making REL_EPS negative effectively shrinks the box
  inline
  auto
  LiangBarskyClip3DBox_plain_types(const T&x_lo, const T&y_lo, const T&z_lo, // input: box lower bounds
				   const T&x_hi, const T&y_hi, const T&z_hi, // input: box upper bounds
				   const T&x_p0, const T&y_p0, const T&z_p0, // input: starting point
				   const T&x_p1, const T&y_p1, const T&z_p1, // input: ending point
				   T& t_min, T&t_max) // output: the line paramter coordinates of the clipping endpoints 
  {
    return LiangBarskyClip3DBox<T,REL_EPS>(Vector3D<T>{{x_lo,y_lo,z_lo}},
					   Vector3D<T>{{x_hi,y_hi,z_hi}},
					   Vector3D<T>{{x_p0,y_p0,z_p0}},
					   Vector3D<T>{{x_p1,y_p1,z_p1}},
					   t_min, t_max);
  }
  
}

#endif
