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


#ifndef KSN_ANALYTICAL_SOLUTIONS_HH
#define KSN_ANALYTICAL_SOLUTIONS_HH

//#define KKC_DEBUG

#include <iostream>
#include <string>
#include <list>
#include <cmath>
#include <vector>
#include <algorithm>
#include <limits>
#include <tuple>

#include "ClipLiangBarsky3D.hh"

/*

This file provides utilities for computing exact solutions of the KSN "type i" (i.e. no scattering) benchmark problems from:

It requires a companion file called ClipLiangBarsky3D.hh.

The easiest way to use this code is something like:
--------------------------------------------------
#include "KSN_AnalyticalSolutions.hh"

using namespace KSN_AnalyticalSolutions;

Problem_3<double> problem;
double x = -4.5;
double y = -10.5;
double z = -9.5;
double xdir = -0.14451;
double ydir = -0.014233;
double zdir = -0.989401;

double flux = compute_exact_solution<double>(x,y,z,
                                             xdir, ydir, zdir,
					     problem);
--------------------------------------------------
where x,y and z are the coordinates where you would like the solution evaluated and 
and xdir, ydir, zdir are components of the unit vector pointing in the streaming direction.

You can also ask for the integrated flux at a point (not normalized by 4pi) by doing something like:
--------------------------------------------------
double integral = integrate_total_flux<double>(x,y,z,
    					       problem,
					       N_azimuthal,N_polar);
--------------------------------------------------
where N_azimuthal and N_polar are the number of quadrature cells in the azimuthal and polar directions.  Note 
4-pt Gauss quadrature is used for the surface integral.
 */


namespace KSN_AnalyticalSolutions {

  enum RegionType {
    Source = 0,
    Void   = 1,
    Shield = 2,
    NumberOfRegionTypes
  };

  // a region keeps track of a 3D box, basically.  A collection of these boxes is used to define each KSN test problem's domain.
  template<typename T>
  struct RegionDefinition {
    RegionDefinition(const T & x_mn, const T &x_mx,
		     const T & y_mn, const T &y_mx,
		     const T & z_mn, const T &z_mx,
		     const RegionType &type,
		     const std::string &desc) : x_min(x_mn), x_max(x_mx),
						y_min(y_mn), y_max(y_mx),
						z_min(z_mn), z_max(z_mx),
						region_type(type),
						description(desc) {}
    
    RegionDefinition(const RegionDefinition<T> &rd) : x_min(rd.x_min), x_max(rd.x_max),
						      y_min(rd.y_min), y_max(rd.y_max),
						      z_min(rd.z_min), z_max(rd.z_max),
						      region_type(rd.region_type),
						      description(rd.description) {}
    
    virtual ~RegionDefinition() {}

  public:
    T x_min, x_max, y_min, y_max, z_min, z_max;
    RegionType region_type;
    std::string description;
    
  private:
    inline RegionDefinition() {};

  };

  template<typename T>
  class TestProblem {
  public:
    // this is basically a class that holds the list of regions for each test problem
    // Note that we define a region as a box, each box can have boxes nexted in it from other regions, the nested regions are effectively carved out of
    //   the enclosing region when the ray-trace is performed.
    // HOWEVER: this is not a generic constructive solid geometry code, it has some shortcuts specific to this application because the developer is lazy
    //          do not think you can simply lay out a bunch of arbitrarily overlapping boxes and expect it to work...
    using iterator = typename std::list<KSN_AnalyticalSolutions::RegionDefinition<T> >::iterator;
    using const_iterator = typename std::list<KSN_AnalyticalSolutions::RegionDefinition<T> >::const_iterator;

    TestProblem() {}
    virtual ~TestProblem() {}
    std::list<KSN_AnalyticalSolutions::RegionDefinition<T> > regions;
  };

  template<typename T>
  class Problem_TEST_GEOM : public TestProblem<T> {
  public:
    Problem_TEST_GEOM()
    {
      this->regions.push_back(RegionDefinition<T>(0, 10,
						  0, 10,
						  0, 10,
						  Source,
						  "Source"));
      
      this->regions.push_back(RegionDefinition<T>(0, 50,
						  0, 50,
						  0, 50,
						  Void,
						  "Void"));
      
      this->regions.push_back(RegionDefinition<T>(0, 100,
						  0, 100,
						  0, 100,
						  Shield,
						  "Shield"));
    }
    
  };


    template<typename T>
  class Problem_TEST_GEOM_2 : public TestProblem<T> {
  public:
    Problem_TEST_GEOM_2()
    {
      this->regions.push_back(RegionDefinition<T>(0, 10,
						  0, 10,
						  0, 10,
						  Source,
						  "Source"));
      
      this->regions.push_back(RegionDefinition<T>(0, 10,
						  10, 50,
						  0, 10,
						  Void,
						  "Void"));
      
      this->regions.push_back(RegionDefinition<T>(0, 10,
						  50, 100,
						  0, 10,
						  Void,
						  "Void"));
      
      this->regions.push_back(RegionDefinition<T>(0, 100,
						  0, 100,
						  0, 100,
						  Shield,
						  "Shield"));
    }
    
  };


  // The test problems define reflection boundary condtions at the x=0, y=0 and z=0 planes.  We just compute it over the whole (unreflected) geometry.
  template<typename T>
  class Problem_1 : public TestProblem<T> {
  public:
    Problem_1()
    {
      this->regions.push_back(RegionDefinition<T>(-10, 10,
						  -10, 10,
						  -10, 10,
						  Source,
						  "Source"));
      
      this->regions.push_back(RegionDefinition<T>(-50, 50,
						  -50, 50,
						  -50, 50,
						  Void,
						  "Void"));
      
      this->regions.push_back(RegionDefinition<T>(-100, 100,
						  -100, 100,
						  -100, 100,
						  Shield,
						  "Shield"));
    }
    
  };

  template<typename T>
  class Problem_2 : public TestProblem<T> {
  public:
    Problem_2()
    {
      this->regions.push_back(RegionDefinition<T>(-10, 10,
						  -10, 10,
						  -10, 10,
						  Source,
						  "Source"));
      
      this->regions.push_back(RegionDefinition<T>(-10, 10,
						  -100,-10,
						  -10, 10,
						  Void,
						  "Void"));

      this->regions.push_back(RegionDefinition<T>(-10, 10,
						   10, 100,
						  -10, 10,
						  Void,
						  "Void"));

      this->regions.push_back(RegionDefinition<T>(-60, 60,
						  -100, 100,
						  -60, 60,
						  Shield,
						  "Shield"));
    }
    
  };

    template<typename T>
  class Problem_3 : public TestProblem<T> {
  public:
    Problem_3()
    {
      this->regions.push_back(RegionDefinition<T>(-10, 10,
						  -10, 10,
						  -10, 10,
						  Source,
						  "Source"));

      // first leg of the channel, along the y axis
      this->regions.push_back(RegionDefinition<T>(-10, 10,
						   10, 50,
						  -10, 10,
						  Void,
						  "Void"));
      this->regions.push_back(RegionDefinition<T>(-10, 10,
						  -50,-10,
						  -10, 10,
						  Void,
						  "Void"));
      // second leg of the channel, parallel to x axis
      this->regions.push_back(RegionDefinition<T>(-40, 40,
						   50, 60,
						  -10, 10,
						  Void,
						  "Void"));

      this->regions.push_back(RegionDefinition<T>(-40, 40,
						  -60,-50,
						  -10, 10,
						  Void,
						  "Void"));

      // third leg of the channel, there is one in each octant
      this->regions.push_back(RegionDefinition<T>(30,40,
						  50,60,
						  10,30,
						  Void,
						  "Void"));

      this->regions.push_back(RegionDefinition<T>(-40,-30,
						   50,60,
						   10,30,
						  Void,
						  "Void"));

      this->regions.push_back(RegionDefinition<T>(30,40,
						  50,60,
						 -30,-10,
						  Void,
						  "Void"));

      this->regions.push_back(RegionDefinition<T>(-40,-30,
						  50,60,
						  -30,-10,
						  Void,
						  "Void"));

      // - negative y legs
      this->regions.push_back(RegionDefinition<T>(30,40,
						  -60,-50,
						  10,30,
						  Void,
						  "Void"));

      this->regions.push_back(RegionDefinition<T>(-40,-30,
						  -60,-50,
						   10,30,
						  Void,
						  "Void"));

      this->regions.push_back(RegionDefinition<T>(30,40,
						  -60,-50,
						  -30,-10,
						  Void,
						  "Void"));

      this->regions.push_back(RegionDefinition<T>(-40,-30,
						  -60,-50,
						  -30,-10,
						  Void,
						  "Void"));
      

      // fourth leg of the channel, also one in each octant
      this->regions.push_back(RegionDefinition<T>(30,40,
						  50,100,
						  30,40,
						  Void,
						  "Void"));
			      
      this->regions.push_back(RegionDefinition<T>(-40,-30,
						  50,100,
						  30,40,
						  Void,
						  "Void"));

      this->regions.push_back(RegionDefinition<T>(30,40,
						  50,100,
						  -40,-30,
						  Void,
						  "Void"));
			      
      this->regions.push_back(RegionDefinition<T>(-40,-30,
						  50,100,
						  -40,-30,
						  Void,
						  "Void"));

      // negative y axis part of fourth channel
      this->regions.push_back(RegionDefinition<T>(30,40,
						  -100,-50,
						  30,40,
						  Void,
						  "Void"));
			      
      this->regions.push_back(RegionDefinition<T>(-40,-30,
						  -100,-50,
						  30,40,
						  Void,
						  "Void"));

      this->regions.push_back(RegionDefinition<T>(30,40,
						  -100,-50,
						  -40,-30,
						  Void,
						  "Void"));
			      
      this->regions.push_back(RegionDefinition<T>(-40,-30,
						  -100,-50,
						  -40,-30,
						  Void,
						  "Void"));

			      
						 
      // Shield, there is only one region
      this->regions.push_back(RegionDefinition<T>(-60, 60,
						  -100, 100,
						  -60, 60,
						  Shield,
						  "Shield"));
    }
    
  };

  template<typename T>
  auto
  compute_exact_solution(const T &x_start, const T &y_start, const T &z_start, //input: x,y,z coordinates to compute the solution at
			 const T &dir_x, const T& dir_y, const T &dir_z, // direction of the flux
			 const TestProblem<T> &test_problem)    //input: the specific test problem to use
  {
    // return output: the angular flux at the given point and direction in units of number/cm^2
    using namespace KSN_AnalyticalSolutions;

    const T source[NumberOfRegionTypes] = {1.0,  0.0, 0.0}; // units are number/(cm^3 s)
    const T sigmat[NumberOfRegionTypes] = {0.1, 1e-4, 0.1}; // units are 1/cm

    // first we loop through the regions to find the longest length scale, we will use this for our starting ray
    T len_scale = 0.0;
    for (typename TestProblem<T>::const_iterator region=test_problem.regions.begin();
	 region!=test_problem.regions.end();
	 region++)
      {
	len_scale = std::max(len_scale,
			    (region->x_max-region->x_min)*(region->x_max-region->x_min) +
			    (region->y_max-region->y_min)*(region->y_max-region->y_min) +
			    (region->z_max-region->z_min)*(region->z_max-region->z_min) );
      }
    len_scale = std::sqrt(len_scale);

    #ifdef KKC_DEBUG
    std::cout<<"max len = "<<len_scale<<std::endl;
    #endif
    
    // compute an endpoint for the ray we are casting using the specified direction
    T x_end = x_start - len_scale * dir_x;//std::cos(theta) * std::cos(phi);
    T y_end = y_start - len_scale * dir_y;//std::sin(theta) * std::cos(phi);
    T z_end = z_start - len_scale * dir_z;//std::sin(phi);

    #ifdef KKC_DEBUG
    std::cout<<"casting ray : ";
    std::cout<<"x("<<x_start<<"->"<<x_end<<")";
    std::cout<<"y("<<y_start<<"->"<<y_end<<")";
    std::cout<<"z("<<z_start<<"->"<<z_end<<")"<<std::endl;
    #endif
    
    struct Segment{ // little helper struct for recording segments and sorting them
      Segment(T &ts, T &te, RegionType rt) : t_s(ts), t_e(te), type(rt) {}
      Segment(const Segment &seg) : t_s(seg.t_s), t_e(seg.t_e), type(seg.type) {}
      T t_s, t_e;
      RegionType type;
    };

    std::vector<Segment> segments;
    
    // loop through all the regions and compute an initial set of segments.  Note these segments will generally overlap, so we fix that up afterwards
    for (typename TestProblem<T>::const_iterator region=test_problem.regions.begin();
	 region!=test_problem.regions.end();
	 region++)
      {
	T t_min=0;
	T t_max=1;
	
	bool segment_inside = KKC_CLIP_3D::LiangBarskyClip3DBox_plain_types<T>(region->x_min, region->y_min, region->z_min,
									       region->x_max, region->y_max, region->z_max,
									       x_start, y_start, z_start,
									       x_end, y_end, z_end,
									       t_min, t_max);
	if (segment_inside)
	  {
	    #ifdef KKC_DEBUG
	    std::cout<<"found a segment with region type "<<region->region_type<<", t_min = "<<t_min<<", t_max = "<<t_max<<std::endl;
	    #endif
	    segments.push_back(Segment(t_min,t_max,region->region_type));
	  }
	
      }

    if (!segments.size())
      return T{0.0};

    std::vector<T> points_on_the_line;
    points_on_the_line.reserve(segments.size()*2);
    for ( auto && segment : segments ) {
      points_on_the_line.push_back(segment.t_s);
      points_on_the_line.push_back(segment.t_e);
    }

    std::sort(points_on_the_line.begin(), points_on_the_line.end()); // I'm not providing a predicate here because I'll take care of precision issues in unique

    auto last_entry = std::unique(points_on_the_line.begin(), points_on_the_line.end(),
				  [](const T&a, const T&b) -> bool
				  {
				    return std::fabs(a-b)<std::numeric_limits<T>::epsilon();// note we use unscaled epsilon because generally a and b are between 0 and 1 and typically around O(0.1) to O(1)
				  }
				  );
    
    points_on_the_line.erase(last_entry, points_on_the_line.end());

    //    for ( T p : points_on_the_line ) std::cout<<p<<std::endl;
    
    std::vector<RegionType> segment_region;
    segment_region.reserve(points_on_the_line.size()-1);
    for ( size_t i=0; i<points_on_the_line.size()-1; i++ )
      {
	const T &mid_pt = 0.5*(points_on_the_line[i] + points_on_the_line[i+1]);
	RegionType region = Shield;
	for ( auto && segment : segments )
	  {
	    if (mid_pt>segment.t_s && mid_pt<segment.t_e)
	      region = std::min(region, segment.type); // there is an assumed priority here, lowest is Shield, highest is Source
	  }
	segment_region.push_back(region);
      }
    
    segments.clear();
    segments.reserve(segment_region.size());
    for ( size_t i=0; i<points_on_the_line.size()-1; i++ )
      {
	segments.push_back(Segment(points_on_the_line[i], points_on_the_line[i+1], segment_region[i]));
      }
    
    #ifdef KKC_DEBUG
    std::cout<<"FINAL SEGMENT LIST"<<std::endl;
    for (typename std::vector<Segment>::iterator i=segments.begin();
	 i!=segments.end();
	 i++)
      {
	std::cout<<"segment t_min="<<i->t_s<<", t_max = "<<i->t_e<<", type = "<<i->type<<std::endl;
      }
    #endif

    if (!segments.size() )
      return T{0.0};
    
    // at this point we have the ray split into segments corresponding to source, void and sheild portions, now its time to do the integral...
    T S = source[segments[0].type];
    T sig = sigmat[segments[0].type];
    T len_seg = len_scale*(segments[0].t_e - segments[0].t_s);
    T integral = (S/sig) * (1.0 - std::exp(-sig*len_seg));
    
    for ( size_t i=1; i<segments.size(); i++ )
      {
	T S_i = source[segments[i].type];
	if ( S_i>std::numeric_limits<T>::epsilon() ) {
	  T sig_i = sigmat[segments[i].type];
	  T sigma_sum = 0.0;
	  for ( size_t k=0; k<i; k++ )
	    {
	      T len_seg = len_scale*(segments[k].t_e - segments[k].t_s);
	      sigma_sum += (sigmat[segments[k].type] - sigmat[segments[i].type]) * len_seg;
	    }
	  integral += (-S_i/sig_i)*std::exp(-sigma_sum)*( std::exp(-sig_i * len_scale * segments[i].t_e) -
							  std::exp(-sig_i * len_scale * segments[i].t_s));
	}
      }
    
    return integral;
  }

  template<typename T>
  inline
  auto
  compute_exact_solution(const T &x_start, const T &y_start, const T &z_start, //input: x,y,z coordinates to compute the solution at
			 const T &theta, const T&phi,        //input: azimuthal and polar coordinates (in radians) for the direction to compute the solution with
			 const TestProblem<T> &test_problem)    //input: the specific test problem to use
  {
    // NOTE: the polar angle, phi, is *measured from the equator* so \phi \in [-pi/2, pi/2]
    //       Why, you ask?  Because.
    const T& dir_x = std::cos(theta) * std::cos(phi);;
    const T& dir_y = std::sin(theta) * std::cos(phi);;
    const T& dir_z = std::sin(phi);

    return compute_exact_solution(x_start, y_start, z_start, dir_x, dir_y, dir_z, test_problem);
  }
  

  template<typename T>
  auto
  integrate_total_flux(const T &x, const T &y, const T &z, //input: x,y,z coordinates to compute the solution at
		       const TestProblem<T> &problem,//input: the specific test problem to use
		       const size_t &N_theta, const size_t &N_phi)    //input: the number of integration intervals in each angular direction
  {
    // return output: an approximation of the total intergated flux

    // we use 4-point gauss quadrature 

    const T pi = std::acos(-1); 
    
    T d_theta = 2*pi/N_theta;
    T d_phi = pi/N_phi;

    const T x1 = std::sqrt(525.0 - 70.0*sqrt(30.))/35.0;
    const T x2 = std::sqrt(525.0 + 70.0*sqrt(30.))/35.0;

    const T w1 = (18.0 + std::sqrt(30))/36.0;
    const T w2 = (18.0 - std::sqrt(30))/36.0;

    constexpr int N = 4;
    T weights[] = { w2,  w1, w1, w2 };
    T abscissa[] = { -x2, -x1, x1, x2 };
    T integral = 0;
    for ( size_t i_theta = 0; i_theta<N_theta; i_theta++)
      for ( size_t i_phi = 0; i_phi<N_phi; i_phi++)
	{
	  T theta[N];
	  T phi[N];
	  for ( size_t i=0; i<N; i++ )
	    {
	      theta[i] = i_theta*d_theta + ((abscissa[i]+1.0)/2.0)*d_theta;
	      phi[i] = -pi/2 + i_phi*d_phi + ((abscissa[i]+1.0)/2.0)*d_phi;
	    }

	  for ( size_t i=0; i<N; i++ )
	    for ( size_t j=0; j<N; j++ )
	      {	
		integral += weights[i]*weights[j]*0.25*d_theta*d_phi * std::cos(phi[j])*compute_exact_solution<T>(x,y,z,
														  theta[i], phi[j],
														  problem);
	      }
	}
    return integral;
  }
}

namespace KSN_Type_ii_reference_solutions {
  // These are NOT exact solutions, but the reported reference results at specific points given in the original paper, computed with a monte carlo method

  template<typename T>
  class TestProblem {
  public:
    using SolutionPoint = std::tuple<T,T,T,T>;  
    using iterator = typename std::vector<SolutionPoint>;
    using const_iterator = typename std::vector<SolutionPoint>;

    TestProblem() {}
    virtual ~TestProblem() {}
    std::vector<SolutionPoint> points;
  };

  template<typename T>
  class Problem_3 : public KSN_Type_ii_reference_solutions::TestProblem<T> {
  public:
    Problem_3() : TestProblem<T>::points{
			 {5, 5,5,8.61578e0},
			 {5,15,5,2.16130e0},
			 {5,25,5,8.93784e-1},
			 {5,35,5,4.78052e-1},
			 {5,45,5,2.89424e-1},
			 {5,55,5,1.92698e-1},
			 {5,65,5,1.04982e-1},
			 {5,75,5,3.37544e-2},
			 {5,85,5,1.08158e-2},
			 {5,95,5,3.39632e-3},
			 
			 { 5,55,5,1.92698e-1},
			 {15,55,5,6.72147e-2},
			 {25,55,5,2.21799e-2},
			 {35,55,5,9.90646e-3},
			 {45,55,5,3.39066e-3},
			 {55,55,5,1.05629e-3},

			 { 5,95,35,3.44804e-4},
			 {15,95,35,2.91825e-4},
			 {25,95,35,2.05793e-4},
			 {35,95,35,2.62086e-4},
			 {45,95,35,1.05367e-4},
			 {55,95,35,4.44962e-5}
    } {}

  };
}

#endif
