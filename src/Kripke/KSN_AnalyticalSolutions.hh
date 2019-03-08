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

#include "ClipLiangBarsky3D.hh"

namespace KSN_AnalyticalSolutions {

  enum RegionType {
    Source = 0,
    Void   = 1,
    Shield = 2,
    NumberOfRegionTypes
  };
    
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
						  -100, 100,
						  -10, 10,
						  Void,
						  "Void"));
#if 0
      this->regions.push_back(RegionDefinition<T>(-10, 10,
						   10, 100,
						  -10, 10,
						  Void,
						  "Void"));
#endif
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
						   -50, 50,
						  -10, 10,
						  Void,
						  "Void"));
      #if 0
      this->regions.push_back(RegionDefinition<T>(-10, 10,
						  -50,-10,
						  -10, 10,
						  Void,
						  "Void"));
      #endif
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
    // return output: the angular flux at the given point and direction
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
    
    // compute an endpoint for the ray we are casting using the angle coordinates
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

    if (!segments.size() )
      return T(0.0);
    

#ifdef KKC_DEBUG 
    std::cout<<"UNPROCESSED INITIAL SEGMENT LIST AFTER SORT "<<std::endl;
    for (typename std::vector<Segment>::iterator i=segments.begin();
	 i!=segments.end();
	 i++)
      {
	std::cout<<"segment t_min="<<i->t_s<<", t_max = "<<i->t_e<<", type = "<<i->type<<std::endl;
      }
#endif
    // ok get an initial sorting of the segments, they might be overlapping so we now have to fix that
    std::sort(segments.begin(), segments.end(), [] (const Segment &a, const Segment &b)->bool
	      {
		if ( std::numeric_limits<T>::epsilon()<std::abs(a.t_s - b.t_s) )
		  return a.t_s<b.t_s;
		else
		  return a.t_e<b.t_e;
	      });

#ifdef KKC_DEBUG 
    std::cout<<"INITIAL SEGMENT LIST AFTER SORT "<<std::endl;
    for (typename std::vector<Segment>::iterator i=segments.begin();
	 i!=segments.end();
	 i++)
      {
	std::cout<<"segment t_min="<<i->t_s<<", t_max = "<<i->t_e<<", type = "<<i->type<<std::endl;
      }
#endif

    // the segments are now sorted first by the starting point and then by the end points
    // it is possible for multiple segments to start on the same point, so fix that now
    for (size_t s=(segments.size()-1); s>0 ; s--)
      {
	  if ( std::numeric_limits<T>::epsilon() > std::abs(segments[s].t_s-segments[s-1].t_s) ) // note the line is parameterized in [0,1]
	    segments[ s ].t_s = segments[s-1].t_e;
      }
    // ok get an initial sorting of the segments, they might be overlapping so we now have to fix that
    std::sort(segments.begin(), segments.end(), [] (const Segment &a, const Segment &b)->bool
	      {
		if ( std::numeric_limits<T>::epsilon()<std::abs(a.t_s - b.t_s) )
		  return a.t_s<b.t_s;
		else
		  return a.t_e<b.t_e;
	      });

    
    #ifdef KKC_DEBUG 
    std::cout<<"INITIAL SEGMENT LIST"<<std::endl;
    for (typename std::vector<Segment>::iterator i=segments.begin();
	 i!=segments.end();
	 i++)
      {
	std::cout<<"segment t_min="<<i->t_s<<", t_max = "<<i->t_e<<", type = "<<i->type<<std::endl;
      }
    #endif
    
    bool done = !segments.size();
    while (!done) {
      bool no_insertion = true;
      done = true;
      for (size_t s=0; s<(segments.size()-1) && no_insertion; s++ )
	{
	  if ( std::numeric_limits<T>::epsilon() < segments[s].t_e-segments[s+1].t_e ) // note the line is parameterized in [0,1]
	    {
	      no_insertion = false;
	      done = false;
	      T te_s = segments[s].t_e;
	      T te_sp1 = segments[s+1].t_e;
	      segments[s].t_e = segments[s+1].t_s;
	      segments.push_back(Segment(te_sp1, te_s, segments[s].type));
	      std::sort(segments.begin(), segments.end(), [] (const Segment &a, const Segment &b)->bool
	      {
		if ( std::numeric_limits<T>::epsilon()<std::abs(a.t_s - b.t_s) )
		  return a.t_s<b.t_s;
		else
		  return a.t_e<b.t_e;
	      });
	    }
	  else if ( std::numeric_limits<T>::epsilon() > std::abs(segments[s].t_e-segments[s+1].t_e) ) // note the line is parameterized in [0,1]
	    {
	      segments[s].t_e = segments[s+1].t_s;
	    }
	}
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
  auto
  compute_exact_solution(const T &x_start, const T &y_start, const T &z_start, //input: x,y,z coordinates to compute the solution at
			 const T &theta, const T&phi,        //input: azimuthal and polar coordinates (in radians) for the direction to compute the solution with
			 const TestProblem<T> &test_problem)    //input: the specific test problem to use
  {
    const T& dir_x = std::cos(theta) * std::cos(phi);;
    const T& dir_y = std::sin(theta) * std::cos(phi);;
    const T& dir_z = std::sin(phi);

    return compute_exact_solution(x_start, y_start, z_start, dir_x, dir_y, dir_z, test_problem);
  }
  

  template<typename T>
  auto
  integrate_total_flux(const T &x, const T &y, const T &z, //input: x,y,z coordinates to compute the solution at
		       const TestProblem<T> &problem,//input: the specific test problem to use
		       const int &N_theta, const int &N_phi)    //input: the number of integration intervals in each angular direction
  {
    // return output: an approximation of the total intergated flux

    // we use 4-point gauss quadrature to reduce the evaluation time (should be 8th order accurate)

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
    for (int i_theta = 0; i_theta<N_theta; i_theta++)
      for (int i_phi = 0; i_phi<N_phi; i_phi++)
	{
	  T theta[N];
	  T phi[N];
	  for ( int i=0; i<N; i++ )
	    {
	      theta[i] = i_theta*d_theta + ((abscissa[i]+1.0)/2.0)*d_theta;
	      phi[i] = -pi/2 + i_phi*d_phi + ((abscissa[i]+1.0)/2.0)*d_phi;
	    }

	  for (int i=0; i<N; i++ )
	    for (int j=0; j<N; j++ )
	      {	
		integral += weights[i]*weights[j]*0.25*d_theta*d_phi * std::cos(phi[j])*compute_exact_solution<T>(x,y,z,
														  theta[i], phi[j],
														  problem);
	      }
	}
    return integral;
  }
}

#endif
