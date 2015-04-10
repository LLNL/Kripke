#include <Kripke/Grid.h>

#include <Kripke/Input_Variables.h>
#include <Kripke/Layout.h>
#include <Kripke/SubTVec.h>
#include <cmath>
#include <sstream>
#include <mpi.h>

#ifdef KRIPKE_USE_SILO
#include <sys/stat.h>
#include <silo.h>
#include <string.h>
#endif

/**
 * Grid_Data constructor
*/
Grid_Data::Grid_Data(Input_Variables *input_vars)
{
  // Create object to describe processor and subdomain layout in space
  // and their adjacencies
  Layout *layout = createLayout(input_vars);

  // create the kernel object based on nesting
  kernel = createKernel(input_vars->nesting, 3);

  // Create quadrature set (for all directions)
  int total_num_directions = input_vars->num_dirsets_per_octant * input_vars->num_dirs_per_dirset;
  InitDirections(this, input_vars);

  num_direction_sets = 8*input_vars->num_dirsets_per_octant;
  num_directions_per_set = input_vars->num_dirs_per_dirset;
  num_group_sets = input_vars->num_groupsets;
  num_groups_per_set = input_vars->num_groups_per_groupset;
  num_zone_sets = 1;
  for(int dim = 0;dim < 3;++ dim){
    num_zone_sets *= input_vars->num_zonesets_dim[dim];
  }

  legendre_order = input_vars->legendre_order;
  total_num_moments = (legendre_order+1)*(legendre_order+1);

  int num_subdomains = num_direction_sets*num_group_sets*num_zone_sets;

  Nesting_Order nest = input_vars->nesting;

  /* Set ncalls */
  niter = input_vars->niter;

  // setup mapping of moments to legendre coefficients
  moment_to_coeff.resize(total_num_moments);
  int nm = 0;
  for(int n = 0;n < legendre_order+1;++ n){
    for(int m = -n;m <= n; ++ m){
      moment_to_coeff[nm] = n;
      ++ nm;
    }
  }

  // setup cross-sections
  int total_num_groups = num_group_sets*num_groups_per_set;
  sigma_tot.resize(total_num_groups, 0.0);

  // Setup scattering transfer matrix for 2 materials
  double sigs_init[3] = {.05, .00005, 0.05};
  sigs.resize(3);
  for(int mat = 0;mat < 3;++ mat){
    // allocate transfer matrix
    sigs[mat] = new SubTVec(kernel->nestingSigs(), total_num_groups, legendre_order+1, total_num_groups);

    // Set to identity for all moments
    sigs[mat]->clear(0.0);
    for(int l = 0;l < legendre_order+1;++ l){
      for(int g = 0;g < total_num_groups;++ g){
        (*sigs[mat])(g, l, g) = sigs_init[mat];
      }
    }
  }

  // just allocate pointer vectors, we will allocate them below
  ell.resize(num_direction_sets, NULL);
  ell_plus.resize(num_direction_sets, NULL);
  phi.resize(num_zone_sets, NULL);
  phi_out.resize(num_zone_sets, NULL);

  // Initialize Subdomains
  zs_to_sdomid.resize(num_zone_sets);
  subdomains.resize(num_subdomains);
  for(int gs = 0;gs < num_group_sets;++ gs){
    for(int ds = 0;ds < num_direction_sets;++ ds){
      for(int zs = 0;zs < num_zone_sets;++ zs){
        // Compupte subdomain id
        int sdom_id = layout->setIdToSubdomainId(gs, ds, zs);

        // Setup the subdomain
        Subdomain &sdom = subdomains[sdom_id];
        sdom.setup(sdom_id, input_vars, gs, ds, zs, directions, kernel, layout);

        // Create ell and ell_plus, if this is the first of this ds
        bool compute_ell = false;
        if(ell[ds] == NULL){
          ell[ds] = new SubTVec(kernel->nestingEll(), total_num_moments, sdom.num_directions, 1);
          ell_plus[ds] = new SubTVec(kernel->nestingEllPlus(), total_num_moments, sdom.num_directions, 1);

          compute_ell = true;
        }

        // Create phi and phi_out, if this is the first of this zs
        if(phi[zs] == NULL){
          zs_to_sdomid[zs] = sdom_id;
          phi[zs] = new SubTVec(nest, total_num_groups, total_num_moments, sdom.num_zones);
          phi_out[zs] = new SubTVec(nest, total_num_groups, total_num_moments, sdom.num_zones);
        }

        // Set the variables for this subdomain
        sdom.setVars(ell[ds], ell_plus[ds], phi[zs], phi_out[zs]);

        if(compute_ell){
          // Compute the L and L+ matrices
          sdom.computeLLPlus();
        }
      }
    }
  }
  delete layout;



  // Now compute number of elements allocated globally
  long long vec_size[4] = {0,0,0,0};
  for(int sdom_id = 0;sdom_id < subdomains.size();++sdom_id){
    Subdomain &sdom = subdomains[sdom_id];
    vec_size[0] += sdom.psi->elements;
    vec_size[1] += sdom.psi->elements;
  }
  for(int zs = 0;zs < num_zone_sets;++ zs){
    vec_size[2] += phi[zs]->elements;
    vec_size[3] += phi_out[zs]->elements;
  }

  long long vec_global[4];
  MPI_Reduce(vec_size, vec_global, 4, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);

  int mpi_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  if(mpi_rank == 0){
    printf("Unknown counts: psi=%ld, rhs=%ld, phi=%ld, phi_out=%ld\n",
      (long)vec_global[0], (long)vec_global[1], (long)vec_global[2], (long)vec_global[3]);
  }
}

Grid_Data::~Grid_Data(){
  delete kernel;
  for(int zs = 0;zs < num_zone_sets;++ zs){
    delete phi[zs];
    delete phi_out[zs];
  }
  for(int ds = 0;ds < num_direction_sets;++ ds){
    delete ell[ds];
    delete ell_plus[ds];
  }
  for(int mat = 0;mat < 3;++ mat){
    delete sigs[mat];
  }
}

/**
 * Randomizes all variables and matrices for testing suite.
 */
void Grid_Data::randomizeData(void){
  for(int i = 0;i < sigma_tot.size();++i){
    sigma_tot[i] = drand48();
  }

  for(int i = 0;i < directions.size();++i){
    directions[i].xcos = drand48();
    directions[i].ycos = drand48();
    directions[i].zcos = drand48();
  }


  for(int s = 0;s < subdomains.size();++ s){
    subdomains[s].randomizeData();
  }

  for(int zs = 0;zs < num_zone_sets;++ zs){
    phi[zs]->randomizeData();
    phi_out[zs]->randomizeData();
  }

  for(int ds = 0;ds < num_direction_sets;++ ds){
    ell[ds]->randomizeData();
    ell_plus[ds]->randomizeData();
  }

  for(int mat = 0;mat < 3;++ mat){
    sigs[mat]->randomizeData();
  }
}

/**
 * Copies all variables and matrices for testing suite.
 * Correctly copies data from one nesting to another.
 */
void Grid_Data::copy(Grid_Data const &b){
  sigma_tot = b.sigma_tot;
  directions = b.directions;

  subdomains.resize(b.subdomains.size());
  for(int s = 0;s < subdomains.size();++ s){
    subdomains[s].copy(b.subdomains[s]);
  }

  for(int zs = 0;zs < num_zone_sets;++ zs){
    phi[zs]->copy(*b.phi[zs]);
    phi_out[zs]->copy(*b.phi_out[zs]);
  }

  for(int ds = 0;ds < ell.size();++ ds){
    ell[ds]->copy(*b.ell[ds]);
    ell_plus[ds]->copy(*b.ell_plus[ds]);
  }

  for(int mat = 0;mat < 3;++ mat){
    sigs[mat]->copy(*b.sigs[mat]);
  }
}

/**
 * Compares all variables and matrices for testing suite.
 * Correctly compares data from one nesting to another.
 */
bool Grid_Data::compare(Grid_Data const &b, double tol, bool verbose){
  bool is_diff = false;

  for(int i = 0;i < directions.size();++i){
    std::stringstream dirname;
    dirname << "directions[" << i << "]";

    is_diff |= compareScalar(dirname.str()+".xcos",
        directions[i].xcos, b.directions[i].xcos, tol, verbose);

    is_diff |= compareScalar(dirname.str()+".ycos",
        directions[i].ycos, b.directions[i].ycos, tol, verbose);

    is_diff |= compareScalar(dirname.str()+".zcos",
        directions[i].zcos, b.directions[i].zcos, tol, verbose);
  }

  for(int s = 0;s < subdomains.size();++ s){
    is_diff |= subdomains[s].compare(
        b.subdomains[s], tol, verbose);

  }
  is_diff |= compareVector("sigma_tot", sigma_tot, b.sigma_tot, tol, verbose);

  for(int zs = 0;zs < num_zone_sets;++ zs){
    is_diff |= phi[zs]->compare("phi", *b.phi[zs], tol, verbose);
    is_diff |= phi_out[zs]->compare("phi_out", *b.phi_out[zs], tol, verbose);
  }

  for(int ds = 0;ds < ell.size();++ ds){
    is_diff |= ell[ds]->compare("ell", *b.ell[ds], tol, verbose);
    is_diff |= ell_plus[ds]->compare("ell_plus", *b.ell_plus[ds], tol, verbose);
  }

  for(int mat = 0;mat < 3;++ mat){
    is_diff |= sigs[mat]->compare("sigs", *b.sigs[mat], tol, verbose);
  }

  return is_diff;
}


#ifdef KRIPKE_USE_SILO
void Grid_Data::writeSilo(std::string const &fname_base){

  // Recompute Phi... so we can write out phi0
  kernel->LTimes(this);


  int mpi_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  int mpi_size;
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  int num_subdomains = subdomains.size();

  if(mpi_rank == 0){
    // Create a root file
    std::string fname_root = fname_base + ".silo";
    DBfile *root = DBCreate(fname_root.c_str(),
        DB_CLOBBER, DB_LOCAL, NULL, DB_HDF5);

    // Write out multimesh for spatial mesh
    {
      // setup mesh names and types
      std::vector<char *> mesh_names(mpi_size*num_zone_sets);
      int mesh_idx = 0;
      for(int rank = 0;rank < mpi_size;++ rank){
        for(int idx = 0;idx < num_zone_sets;++ idx){
          int sdom_id = zs_to_sdomid[idx];
          std::stringstream name;

          name << fname_base << "/rank_" << rank << ".silo:/sdom" << sdom_id << "/mesh";
          mesh_names[mesh_idx] = strdup(name.str().c_str());

          mesh_idx ++;
        }
      }
      std::vector<int> mesh_types(mpi_size*num_zone_sets, DB_QUAD_RECT);

      DBPutMultimesh(root, "mesh", mpi_size*num_zone_sets,
          &mesh_names[0], &mesh_types[0], NULL);

      // cleanup
      for(int i = 0;i < mpi_size*num_zone_sets; ++i){
        free(mesh_names[i]);
      }
    }

    // Write out multimat for materials
    {
      // setup mesh names and types
      std::vector<char *> mat_names(mpi_size*num_zone_sets);
      int mesh_idx = 0;
      for(int rank = 0;rank < mpi_size;++ rank){
        for(int idx = 0;idx < num_zone_sets;++ idx){
          int sdom_id = zs_to_sdomid[idx];
          std::stringstream name;

          name << fname_base << "/rank_" << rank << ".silo:/sdom" << sdom_id << "/material";
          mat_names[mesh_idx] = strdup(name.str().c_str());

          mesh_idx ++;
        }
      }

      DBPutMultimat(root, "material", mpi_size*num_zone_sets,
          &mat_names[0],  NULL);

      // cleanup
      for(int i = 0;i < mpi_size*num_zone_sets; ++i){
        free(mat_names[i]);
      }
    }

    // Write out multivar for phi0
    {
      // setup mesh names and types
      std::vector<char *> var_names(mpi_size*num_zone_sets);
      int mesh_idx = 0;
      for(int rank = 0;rank < mpi_size;++ rank){
        for(int idx = 0;idx < num_zone_sets;++ idx){
          int sdom_id = zs_to_sdomid[idx];
          std::stringstream name;

          name << fname_base << "/rank_" << rank << ".silo:/sdom" << sdom_id << "/phi0";
          var_names[mesh_idx] = strdup(name.str().c_str());

          mesh_idx ++;
        }
      }

      std::vector<int> var_types(mpi_size*num_zone_sets, DB_QUADVAR);

      DBPutMultivar(root, "phi0", mpi_size*num_zone_sets,
          &var_names[0],  &var_types[0] , NULL);

      // cleanup
      for(int i = 0;i < mpi_size*num_zone_sets; ++i){
        free(var_names[i]);
      }
    }

    // Root file
    DBClose(root);

    // Create a subdirectory to hold processor info
    mkdir(fname_base.c_str(), 0750);
  }

  // Sync up, so everyone sees the subdirectory
  MPI_Barrier(MPI_COMM_WORLD);

  // Create our processor file
  std::stringstream ss_proc;
  ss_proc << fname_base << "/rank_" << mpi_rank << ".silo";
  DBfile *proc = DBCreate(ss_proc.str().c_str(),
      DB_CLOBBER, DB_LOCAL, NULL, DB_HDF5);

  // Write out data for each subdomain
  for(int sdom_id = 0;sdom_id < num_subdomains;++ sdom_id){
    Subdomain &sdom = subdomains[sdom_id];

    // Create a directory for the subdomain
    std::stringstream dirname;
    dirname << "/sdom" << sdom_id;
    DBMkDir(proc, dirname.str().c_str());

    // Set working directory
    DBSetDir(proc, dirname.str().c_str());

    // Write the mesh
    {
      char *coordnames[3] = {"X", "Y", "Z"};
      double *coords[3];
      for(int dim = 0;dim < 3;++ dim){
        coords[dim] = new double[sdom.nzones[dim]+1];
        coords[dim][0] = sdom.zeros[dim];
        for(int z = 0;z < sdom.nzones[dim];++ z){
          coords[dim][1+z] = coords[dim][z] + sdom.deltas[dim][z+1];
        }
      }
      int nnodes[3] = {
          sdom.nzones[0]+1,
          sdom.nzones[1]+1,
          sdom.nzones[2]+1
      };

      DBPutQuadmesh(proc, "mesh", coordnames, coords, nnodes, 3, DB_DOUBLE,
          DB_COLLINEAR, NULL);

      // cleanup
      delete[] coords[0];
      delete[] coords[1];
      delete[] coords[2];
    }

    // Write the material
    {
      int num_zones = sdom.num_zones;
      int num_mixed = sdom.mixed_material.size();
      int matnos[3] = {1, 2, 3};
      std::vector<int> matlist(num_zones, 0);
      std::vector<int> mix_next(num_mixed, 0);
      std::vector<int> mix_mat(num_mixed, 0);

      // setup matlist and mix_next arrays
      int last_z = -1;
      for(int m = 0;m < num_mixed;++ m){
        mix_mat[m] = sdom.mixed_material[m] + 1;
        int z = sdom.mixed_to_zones[m];
        if(matlist[z] == 0){
            matlist[z] = -(1+m);
        }
        // if we are still on the same zone, make sure the last mix points
        // here
        if(z == last_z){
          mix_next[m-1] = m+1;
        }
        last_z = z;
      }

      DBPutMaterial(proc, "material", "mesh", 3, matnos,
          &matlist[0], sdom.nzones, 3,
          &mix_next[0], &mix_mat[0], &sdom.mixed_to_zones[0], &sdom.mixed_fraction[0], num_mixed,
          DB_DOUBLE, NULL);
    }

    // Write phi0
    {
      int num_zones = sdom.num_zones;
      std::vector<double> phi0(num_zones);

      // extract phi0 from phi for the 0th group
      for(int z = 0;z < num_zones;++ z){
        phi0[z] = (*sdom.phi)(0,0,z);
      }

      DBPutQuadvar1(proc, "phi0", "mesh", &phi0[0],
          sdom.nzones, 3, NULL, 0, DB_DOUBLE, DB_ZONECENT, NULL);
    }
  }

  // Close processor file
  DBClose(proc);
}
#endif


