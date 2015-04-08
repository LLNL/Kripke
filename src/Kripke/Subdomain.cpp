#include <Kripke/Grid.h>
#include <Kripke/SubTVec.h>
#include <Kripke/Input_Variables.h>

#include <cmath>
#include <sstream>


namespace {
  /**
    This function defined the material distribution in space.
    This defines Problem 3 from Kobayashi
    Where Region 1 is material 0, 2 is 1 and 3 is 2.
  */
  inline int queryMaterial(double x, double y, double z){
    // Problem is defined for one octant, with reflecting boundaries
    // We "unreflect" it here by taking abs values
    x = std::abs(x);
    y = std::abs(y);
    z = std::abs(z);

    // Central 20x20x20 box is Region 1
    if(x <= 10.0 && y <= 10.0 && z <= 10.0){
      return 0;
    }

    // Leg 1 of Region 2
    if(x <= 10.0 && y <= 60.0 && z <= 10.0){
      return 1;
    }

    // Leg 2 of Region 2
    if(x <= 40.0 && y >= 50.0 && y <= 60.0 && z <= 10.0){
      return 1;
    }

    // Leg 3 of Region 2
    if(x >= 30.0 && x <= 40.0 && y >= 50.0 && y <= 60.0 && z <= 40.0){
      return 1;
    }

    // Leg 4 of Region 2
    if(x >= 30.0 && x <= 40.0 && y >= 50.0 && z >= 30.0 && z <= 40.0){
      return 1;
    }

    // Rest is filled with region 3
    return 2;
  }
}



Subdomain::Subdomain() :
  idx_dir_set(0),
  idx_group_set(0),
  idx_zone_set(0),
  num_groups(0),
  num_directions(0),
  num_zones(0),
  group0(0),
  direction0(0),
  psi(NULL),
  rhs(NULL),
  sigt(NULL),
  directions(NULL),
  ell(NULL),
  ell_plus(NULL),
  phi(NULL),
  phi_out(NULL)
{
  plane_data[0] = NULL;
  plane_data[1] = NULL;
  plane_data[2] = NULL;
}
Subdomain::~Subdomain(){
  delete psi;
  delete rhs;
  delete sigt;
  delete plane_data[0];
  delete plane_data[1];
  delete plane_data[2];
}


/**
  Setup subdomain and allocate data
*/
void Subdomain::setup(int sdom_id, Input_Variables *input_vars, int gs, int ds, int zs,
    std::vector<Directions> &direction_list, Kernel *kernel, Layout *layout)
{
  // set the set indices
  idx_group_set = gs;
  idx_dir_set = ds;
  idx_zone_set = zs;

  num_groups = input_vars->num_groups_per_groupset;
  group0 = gs * input_vars->num_groups_per_groupset;

  num_directions = input_vars->num_dirs_per_dirset;
  direction0 = ds * input_vars->num_dirs_per_dirset;
  directions = &direction_list[direction0];

  num_zones = 1;
  for(int dim = 0;dim < 3;++ dim){
    // Compute number of zones in this dimension
    nzones[dim] = layout->getNumZones(sdom_id, dim);
    num_zones *= nzones[dim];

    // Compute grid deltas in this dimension (including ghost zone deltas)
    std::pair<double, double> dim_extent = layout->getSpatialExtents(sdom_id, dim);
    zeros[dim] = dim_extent.first;
    double dx = (dim_extent.second-dim_extent.first)/(double)nzones[dim];
    deltas[dim].resize(nzones[dim]+2);
    for(int z = 0;z < nzones[dim]+2;++ z){
      deltas[dim][z] = dx;
    }
  }

  // allocate storage for the sweep boundary data (upwind and downwind share same buffer)
  plane_data[0] = new SubTVec(kernel->nestingPsi(), num_groups, num_directions, nzones[1] * nzones[2]);
  plane_data[1] = new SubTVec(kernel->nestingPsi(), num_groups, num_directions, nzones[0] * nzones[2]);
  plane_data[2] = new SubTVec(kernel->nestingPsi(), num_groups, num_directions, nzones[0] * nzones[1]);

  // allocate the storage for solution and source terms
  psi = new SubTVec(kernel->nestingPsi(), num_groups, num_directions, num_zones);
  rhs = new SubTVec(kernel->nestingPsi(), num_groups, num_directions, num_zones);
  sigt = new SubTVec(kernel->nestingSigt(), num_groups, 1, num_zones);

  computeSweepIndexSet();

  // Setup neighbor data
  int dirs[3] = { directions[0].id, directions[0].jd, directions[0].kd};
  for(int dim = 0;dim < 3;++ dim){
    downwind[dim] = layout->getNeighbor(sdom_id, dim, dirs[dim]);
    upwind[dim] = layout->getNeighbor(sdom_id, dim, -1 * dirs[dim]);
  }

  // paint the mesh
  int num_subsamples = 4; // number of subsamples per spatial dimension
  double sample_vol_frac = 1.0 / (double)(num_subsamples*num_subsamples*num_subsamples);
  int zone_id = 0;
  double pz = zeros[2];

  for (int k = 0; k < nzones[2]; k++) {
    double sdz = deltas[2][k+1] / (double)(num_subsamples+1);
    double py = zeros[1];

    for (int j = 0; j != nzones[1]; j ++) {
      double sdy = deltas[1][j+1] / (double)(num_subsamples+1);
      double px = zeros[0];

      for (int i = 0; i != nzones[0]; i ++) {
        double sdx = deltas[0][i+1] / (double)(num_subsamples+1);

        // subsample probe the geometry to get our materials
        double frac[3] = {0.0, 0.0, 0.0}; // fraction of both materials
        double spz = pz + sdz;

        for(int sk = 0;sk < num_subsamples;++ sk){
          double spy = py + sdy;
          for(int sj = 0;sj < num_subsamples;++ sj){
            double spx = px + sdx;
            for(int si = 0;si < num_subsamples;++ si){

              int mat = queryMaterial(spx, spy, spz);
              frac[mat] += sample_vol_frac;

              spx += sdx;
            }
            spy += sdy;
          }
          spz += sdz;
        }

        // Add material to zone
        for(int mat = 0;mat < 3;++ mat){
          if(frac[mat] > 0.0){
            mixed_to_zones.push_back(zone_id);
            mixed_material.push_back(mat);
            mixed_fraction.push_back(frac[mat]);
          }
        }

        // increment zone
        px += deltas[0][i+1];
        zone_id ++;
      }
      py += deltas[1][j+1];
    }
    pz += deltas[2][k+1];
  }
}

void Subdomain::setVars(SubTVec *ell_ptr, SubTVec *ell_plus_ptr,
    SubTVec *phi_ptr, SubTVec *phi_out_ptr){

  ell = ell_ptr;
  ell_plus = ell_plus_ptr;
  phi = phi_ptr;
  phi_out = phi_out_ptr;
}


/**
 * Randomizes data for a set.
 */
void Subdomain::randomizeData(void){
  psi->randomizeData();
  rhs->randomizeData();
  sigt->randomizeData();

  for(int d = 0;d < 3;++ d){
    for(int i = 0;i < deltas[d].size();++ i){
      deltas[d][i] = drand48();
    }
  }
}

/**
 * Copies two sets, allowing for different nestings.
 */
void Subdomain::copy(Subdomain const &b){
  psi->copy(*b.psi);
  rhs->copy(*b.rhs);
  sigt->copy(*b.sigt);

  for(int d = 0;d < 3;++ d){
    deltas[d] = b.deltas[d];
  }
}

/**
 * Compares two sets, allowing for different nestings.
 */
bool Subdomain::compare(Subdomain const &b, double tol, bool verbose){
  std::stringstream namess;
  namess << "gdset[gs=" << idx_group_set << ", ds=" << idx_dir_set << ", zs=" << idx_zone_set << "]";
  std::string name = namess.str();

  bool is_diff = false;
  is_diff |= psi->compare(name+".psi", *b.psi, tol, verbose);
  is_diff |= rhs->compare(name+".rhs", *b.rhs, tol, verbose);
  is_diff |= sigt->compare(name+".sigt", *b.sigt, tol, verbose);

  is_diff |= compareVector(name+".deltas[0]", deltas[0], b.deltas[0], tol, verbose);
  is_diff |= compareVector(name+".deltas[1]", deltas[1], b.deltas[1], tol, verbose);
  is_diff |= compareVector(name+".deltas[2]", deltas[2], b.deltas[2], tol, verbose);

  return is_diff;
}

/**
 * Compute sweep index sets.
 * Determines logical indices, and increments for i,j,k based on grid
 * information and quadrature set sweeping direction.
 */
void Subdomain::computeSweepIndexSet(void){
  if(directions[0].id > 0){
    sweep_block.start_i = 0;
    sweep_block.end_i = nzones[0]-1;
    sweep_block.inc_i = 1;
  }
  else {
    sweep_block.start_i = nzones[0]-1;
    sweep_block.end_i = 0;
    sweep_block.inc_i = -1;
  }

  if(directions[0].jd > 0){
    sweep_block.start_j = 0;
    sweep_block.end_j = nzones[1]-1;
    sweep_block.inc_j = 1;
  }
  else {
    sweep_block.start_j = nzones[1]-1;
    sweep_block.end_j = 0;
    sweep_block.inc_j = -1;
  }

  if(directions[0].kd > 0){
    sweep_block.start_k = 0;
    sweep_block.end_k = nzones[2]-1;
    sweep_block.inc_k =  1;
  }
  else {
    sweep_block.start_k = nzones[2]-1;
    sweep_block.end_k = 0;
    sweep_block.inc_k = -1;
  }
  
  
//LG
    sweep_block.Nhyperplanes = 0;
    sweep_block.ii_jj_kk_z_idx = NULL;
    {

      int N = nzones[0];
      if (nzones[1] > N) N=nzones[1];
      if (nzones[2] > N) N=nzones[2];

      int i_inc = sweep_block.inc_i;
      int j_inc = sweep_block.inc_j;
      int k_inc = sweep_block.inc_k;
      int i_min, i_max, j_min, j_max, k_min, k_max;
      int counter = 0;
      int Nslices = 0;
      int offset[3*N+4];
      offset[0] = 0;

      if ( i_inc == 1){
        i_min = sweep_block.start_i;
        i_max = sweep_block.end_i-1;
      }
      else{
        i_min = sweep_block.end_i + 1;
        i_max = sweep_block.start_i;
      }
      if ( j_inc == 1){
        j_min = sweep_block.start_j;
        j_max = sweep_block.end_j-1;
      }
      else{
        j_min = sweep_block.end_j + 1;
        j_max = sweep_block.start_j;
      }
      if ( k_inc == 1){
        k_min = sweep_block.start_k;
        k_max = sweep_block.end_k-1;
      }
      else{
        k_min = sweep_block.end_k + 1;
        k_max = sweep_block.start_k;
      }
      int ii_tmp = (1 - i_inc)/2*i_max;
      int jj_tmp = (1 - j_inc)/2*j_max;
      int kk_tmp = (1 - k_inc)/2*k_max;


      sweep_block.ii_jj_kk_z_idx = new int[nzones[0]*nzones[1]*nzones[2]*4];
      sweep_block.d_ii_jj_kk_z_idx = (int*) get_cudaMalloc(size_t   (nzones[0]*nzones[1]*nzones[2]*4) * sizeof(int)   );

      int *ii_jj_kk_z_idx = sweep_block.ii_jj_kk_z_idx;

      for (int C = 0; C <=(3*N); ++C){   //for each C we can touch zone["i","j","k"]  as well as "d" and "group"    in parallel
       int FLAG=0;
       for (int i = 0; i <= C; ++i){
         for (int j = 0; j <= C; ++j){
            int k = C - i - j; // surface equation i+j+j=C
            //flip if needed

            int ii = ii_tmp + i*i_inc;
            int jj = jj_tmp + j*j_inc;
            int kk = kk_tmp + k*k_inc;

            if (ii <= i_max && jj <= j_max && kk <= k_max && ii >= i_min && jj >= j_min && kk >= k_min){
              ii_jj_kk_z_idx[counter*4] = ii;
              ii_jj_kk_z_idx[counter*4+1] = jj;
              ii_jj_kk_z_idx[counter*4+2] = kk;
              ii_jj_kk_z_idx[counter*4+3] =   ii + nzones[0]*jj + nzones[0]*nzones[1]*kk;//  Zonal_INDEX(ii, jj, kk);
              counter++; //counts all elements
              FLAG++;   //counts elements per slice
           }
         }
       }
       if (FLAG){
          Nslices++;
          offset[Nslices] = offset[Nslices-1] + FLAG;
       }
    }
    sweep_block.Nhyperplanes = Nslices;
    sweep_block.offset = new int[Nslices+1];
    for (int ii = 0; ii <= Nslices; ++ii) sweep_block.offset[ii] = offset[ii];

    do_cudaMemcpyH2D( (void *) (sweep_block.d_ii_jj_kk_z_idx), (void *) ii_jj_kk_z_idx, (size_t) (nzones[0]*nzones[1]*nzones[2]*4) * sizeof(int));
    sweep_block.d_offset = (int*) get_cudaMalloc(size_t (Nslices+1) * sizeof(int));
    do_cudaMemcpyH2D( (void *) (sweep_block.d_offset), (void *) offset, size_t (Nslices+1)  * sizeof(int));

   }
  }
}



