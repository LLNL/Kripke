#include<Kripke/Kernel/Kernel_3d_ZDG.h>
#include<Kripke/Grid.h>
#include<Kripke/SubTVec.h>

Kernel_3d_ZDG::Kernel_3d_ZDG() {

}

Kernel_3d_ZDG::~Kernel_3d_ZDG() {

}

Nesting_Order Kernel_3d_ZDG::nestingPsi(void) const {
  return NEST_ZDG;
}

Nesting_Order Kernel_3d_ZDG::nestingPhi(void) const {
  return NEST_ZDG;
}

Nesting_Order Kernel_3d_ZDG::nestingSigt(void) const {
  return NEST_DZG;
}

Nesting_Order Kernel_3d_ZDG::nestingEll(void) const {
  return NEST_ZDG;
}

Nesting_Order Kernel_3d_ZDG::nestingEllPlus(void) const {
  return NEST_ZDG;
}

Nesting_Order Kernel_3d_ZDG::nestingSigs(void) const {
  return NEST_ZDG;
}


void Kernel_3d_ZDG::LTimes(Grid_Data *grid_data) {
  // Outer parameters
  int nidx = grid_data->total_num_moments;

  // Clear phi
  for(int ds = 0;ds < grid_data->num_zone_sets;++ ds){
    grid_data->phi[ds]->clear(0.0);
  }

  // Loop over Subdomains
  int num_subdomains = grid_data->subdomains.size();
  for (int sdom_id = 0; sdom_id < num_subdomains; ++ sdom_id){
    Subdomain &sdom = grid_data->subdomains[sdom_id];

    // Get dimensioning
    int num_groups = sdom.phi->groups;
    int num_zones = sdom.num_zones;
    int num_local_groups = sdom.num_groups;
    int group0 = sdom.group0;
    int num_local_directions = sdom.num_directions;

    /* 3D Cartesian Geometry */
    double * KRESTRICT ell_d_ptr = sdom.ell->ptr();

#ifdef KRIPKE_USE_OPENMP
#pragma omp parallel for
#endif
    for (int z = 0; z < num_zones; z++) {
      double * KRESTRICT psi = sdom.psi->ptr(0, 0, z);
      double * KRESTRICT ell_d = ell_d_ptr;

      for (int d = 0; d < num_local_directions; d++) {
        double * KRESTRICT phi = sdom.phi->ptr(group0, 0, z);

        for(int nm_offset = 0;nm_offset < nidx;++nm_offset){
          double ell_d_nm = ell_d[nm_offset];

          for (int group = 0; group < num_local_groups; ++group) {
            phi[group] += ell_d_nm * psi[group];
          }
          phi += num_groups;
        }
        ell_d += nidx;
        psi += num_local_groups;
      }
    }

  } // Subdomain
}

void Kernel_3d_ZDG::LPlusTimes(Grid_Data *grid_data) {
  // Outer parameters
  int nidx = grid_data->total_num_moments;

  // Loop over Subdomains
  int num_subdomains = grid_data->subdomains.size();
  for (int sdom_id = 0; sdom_id < num_subdomains; ++ sdom_id){
    Subdomain &sdom = grid_data->subdomains[sdom_id];

    // Get dimensioning
    int num_zones = sdom.num_zones;
    int num_groups = sdom.phi_out->groups;
    int num_local_groups = sdom.num_groups;
    int group0 = sdom.group0;
    int num_local_directions = sdom.num_directions;

    // Get Variables
    sdom.rhs->clear(0.0);

    /* 3D Cartesian Geometry */
    double * KRESTRICT ell_plus_ptr = sdom.ell_plus->ptr();

#ifdef KRIPKE_USE_OPENMP
#pragma omp parallel for
#endif
    for (int z = 0; z < num_zones; z++) {
      double * KRESTRICT rhs = sdom.rhs->ptr(0, 0, z);

      double * ell_plus_d = ell_plus_ptr;
      for (int d = 0; d < num_local_directions; d++) {

        double * KRESTRICT phi_out = sdom.phi_out->ptr(group0, 0, z);

        for(int nm_offset = 0;nm_offset < nidx;++nm_offset){
          double ell_plus_d_n_m = ell_plus_d[nm_offset];

          for (int group = 0; group < num_local_groups; ++group) {
            rhs[group] += ell_plus_d_n_m * phi_out[group];
          }
          phi_out += num_groups;
        }
        rhs += num_local_groups;
        ell_plus_d += nidx;
      }
    }
  } // Subdomain
}


/**
  Compute scattering source term phi_out from flux moments in phi.
*/
void Kernel_3d_ZDG::scattering(Grid_Data *grid_data){

}

/* Sweep routine for Diamond-Difference */
/* Macros for offsets with fluxes on cell faces */
#define I_PLANE_INDEX(j, k) (k)*(local_jmax) + (j)
#define J_PLANE_INDEX(i, k) (k)*(local_imax) + (i)
#define K_PLANE_INDEX(i, j) (j)*(local_imax) + (i)
#define Zonal_INDEX(i, j, k) (i) + (local_imax)*(j) \
  + (local_imax)*(local_jmax)*(k)

void Kernel_3d_ZDG::sweep(Subdomain *sdom) {
  int num_directions = sdom->num_directions;
  int num_groups = sdom->num_groups;
  int num_zones = sdom->num_zones;

  Directions *direction = sdom->directions;

  int local_imax = sdom->nzones[0];
  int local_jmax = sdom->nzones[1];
  int local_kmax = sdom->nzones[2];

  double * dx = &sdom->deltas[0][0];
  double * dy = &sdom->deltas[1][0];
  double * dz = &sdom->deltas[2][0];

  // Upwind/Downwind face flux data
  SubTVec &i_plane = *sdom->plane_data[0];
  SubTVec &j_plane = *sdom->plane_data[1];
  SubTVec &k_plane = *sdom->plane_data[2];

  // All directions have same id,jd,kd, since these are all one Direction Set
  // So pull that information out now
  Grid_Sweep_Block const &extent = sdom->sweep_block;

  for (int k = extent.start_k; k != extent.end_k; k += extent.inc_k) {
    double dzk = dz[k + 1];
    for (int j = extent.start_j; j != extent.end_j; j += extent.inc_j) {
      double dyj = dy[j + 1];
      for (int i = extent.start_i; i != extent.end_i; i += extent.inc_i) {
        double dxi = dx[i + 1];

        int z = Zonal_INDEX(i, j, k);
        double * KRESTRICT sigt_z = sdom->sigt->ptr(0, 0, z);

#ifdef KRIPKE_USE_OPENMP
#pragma omp parallel for
#endif
        for (int d = 0; d < num_directions; ++d) {
          double xcos = direction[d].xcos;
          double ycos = direction[d].ycos;
          double zcos = direction[d].zcos;

          double zcos_dzk = 2.0 * zcos / dzk;
          double ycos_dyj = 2.0 * ycos / dyj;
          double xcos_dxi = 2.0 * xcos / dxi;

          double * KRESTRICT psi_z_d = sdom->psi->ptr(0, d, z);
          double * KRESTRICT rhs_z_d = sdom->rhs->ptr(0, d, z);

          double * KRESTRICT psi_lf_z_d = i_plane.ptr(0, d, I_PLANE_INDEX(j, k));
          double * KRESTRICT psi_fr_z_d = j_plane.ptr(0, d, J_PLANE_INDEX(i, k));
          double * KRESTRICT psi_bo_z_d = k_plane.ptr(0, d, K_PLANE_INDEX(i, j));

          for (int group = 0; group < num_groups; ++group) {
            /* Calculate new zonal flux */
            double psi_z_d_g = (rhs_z_d[group]
                + psi_lf_z_d[group] * xcos_dxi
                + psi_fr_z_d[group] * ycos_dyj
                + psi_bo_z_d[group] * zcos_dzk)
                / (xcos_dxi + ycos_dyj + zcos_dzk + sigt_z[group]);

            psi_z_d[group] = psi_z_d_g;

            /* Apply diamond-difference relationships */
            psi_z_d_g *= 2.0;
            psi_lf_z_d[group] = psi_z_d_g - psi_lf_z_d[group];
            psi_fr_z_d[group] = psi_z_d_g - psi_fr_z_d[group];
            psi_bo_z_d[group] = psi_z_d_g - psi_bo_z_d[group];
          }
        }
      }
    }
  }
}

