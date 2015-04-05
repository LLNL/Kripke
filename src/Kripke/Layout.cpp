#include<Kripke/Layout.h>

#include<Kripke/Input_Variables.h>
#include<mpi.h>

Layout::Layout(Input_Variables *input_vars){
  num_group_sets = input_vars->num_groupsets;
  num_direction_sets = input_vars->num_dirsets_per_octant * 8;
  num_zone_sets = 1;
  num_zone_sets_dim[0] = 1;
  num_zone_sets_dim[1] = 1;
  num_zone_sets_dim[2] = 1;

  // grab total number of zones
  total_zones[0] = input_vars->nx;
  total_zones[1] = input_vars->ny;
  total_zones[2] = input_vars->nz;

  // Grab size of processor grid
  num_procs[0] = input_vars->npx;
  num_procs[1] = input_vars->npy;
  num_procs[2] = input_vars->npz;

  /* Set the requested processor grid size */
  int R = num_procs[0] * num_procs[1] * num_procs[2];

  /* Check requested size is the same as MPI_COMM_WORLD */
  int size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if(R != size){
    int myid;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    if(myid == 0){
      printf("ERROR: Incorrect number of MPI tasks. Need %d MPI tasks.", R);
    }
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  /* Compute the local coordinates in the processor decomposition */
  int mpi_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  our_rank[0] = mpi_rank % num_procs[0];
  our_rank[1] = ((mpi_rank - our_rank[0]) / num_procs[0]) % num_procs[1];
  our_rank[2] = (mpi_rank - our_rank[0] - num_procs[0]*our_rank[1]) / (num_procs[0] * num_procs[1]);

  /* Compute the processor neighbor array given the above ordering */
  // mynbr[dimension][direction] = mpi rank of neighbor
  // dimension: 0=x, 1=y, 2=z
  // direction: 0=positive, 1=negative
  // a rank of -1 is used to denote a boundary contition
  mynbr[0][0] = (our_rank[0] == 0)              ? -1 : mpi_rank-1;
  mynbr[0][1] = (our_rank[0] == num_procs[0]-1) ? -1 : mpi_rank+1;
  mynbr[1][0] = (our_rank[1] == 0)              ? -1 : mpi_rank-num_procs[0];
  mynbr[1][1] = (our_rank[1] == num_procs[1]-1) ? -1 : mpi_rank+num_procs[0];
  mynbr[2][0] = (our_rank[2] == 0)              ? -1 : mpi_rank-num_procs[0]*num_procs[1];
  mynbr[2][1] = (our_rank[2] == num_procs[2]-1) ? -1 : mpi_rank+num_procs[0]*num_procs[1];
}
Layout::~Layout(){

}

/**
  Computes the subdomain ID based on a given groupset, directionset, and zoneset.
*/
int Layout::setIdToSubdomainId(int gs, int ds, int zs) const{
  int sdom_id = gs*num_direction_sets*num_zone_sets +
                ds*num_zone_sets +
                zs;

  return sdom_id;
}

/**
  Computes groupset, directionset, and zoneset from a subdomain ID.
*/
void Layout::subdomainIdToSetId(int sdom_id, int &gs, int &ds, int &zs) const {
  gs = sdom_id / (num_direction_sets*num_zone_sets);
  sdom_id = sdom_id % (num_direction_sets*num_zone_sets);

  ds = sdom_id / num_zone_sets;

  zs = sdom_id % num_zone_sets;
}

/**
  Computes the zoneset id along a particular dimension.
*/
int Layout::subdomainIdToZoneSetDim(int sdom_id, int dim) const{
  // Compute zoneset
  int gs, ds, zs;
  subdomainIdToSetId(sdom_id, gs, ds, zs);

  // Compute zone set
  int zs_dim[3];
  zs_dim[0] = zs % num_zone_sets_dim[0];
  zs = zs / num_zone_sets_dim[0];
  zs_dim[1] = zs % num_zone_sets_dim[1];
  zs_dim[2] = zs / num_zone_sets_dim[1];

  return zs_dim[dim];
}


BlockLayout::BlockLayout(Input_Variables *input_vars) :
  Layout(input_vars)
{

}
BlockLayout::~BlockLayout(){

}

Neighbor BlockLayout::getNeighbor(int our_sdom_id, int dim, int dir) const{
  Neighbor n;
  n.mpi_rank = mynbr[dim][dir];
  n.subdomain_id = our_sdom_id;

  return n;
}

/**
  Compute the spatial extents of a subdomain along a given dimension.
*/
std::pair<double, double> BlockLayout::getSpatialExtents(int sdom_id, int dim) const{

  // Start with global problem dimensions
  std::pair<double, double> ext_global(-1.0, 1.0);

  // Subdivide by number of processors in specified dimension
  double dx = (ext_global.second - ext_global.first) / (double)num_procs[dim];
  std::pair<double, double> ext_proc(
    ext_global.first + dx*(double)our_rank[dim],
    ext_global.first + dx*(double)(our_rank[dim] + 1)
  );

  // get the zoneset index along the specified dimension
  int zs_dim = subdomainIdToZoneSetDim(sdom_id, dim);

  // Subdivide by number of subdomains in specified dimension
  double sdx = (ext_proc.second - ext_proc.first) / (double)num_zone_sets_dim[dim];
  std::pair<double, double> ext_sdom(
    ext_proc.first + sdx*(double)zs_dim,
    ext_proc.first + sdx*(double)(zs_dim + 1)
  );

  return ext_sdom;
}

/**
  Computes the number of zones in this subdomain, along specified dimension.
*/
int BlockLayout::getNumZones(int sdom_id, int dim) const{

  // get the zoneset index along the specified dimension
  int zs_dim = subdomainIdToZoneSetDim(sdom_id, dim);

  int total_subdomains = num_procs[dim] * num_zone_sets_dim[dim];
  int global_subdomain  = num_zone_sets_dim[dim] * our_rank[dim] + zs_dim;

  // Compute subset of global zone indices
  int num_zones = total_zones[dim] / total_subdomains;
  int rem = total_zones[dim] % total_subdomains;
  if(rem != 0 && global_subdomain < rem){
    num_zones ++;
  }

  return num_zones;
}

/**
  Factory to create Layout object based on user defined inputs
*/
Layout *createLayout(Input_Variables *input_vars){
  return new BlockLayout(input_vars);
}
