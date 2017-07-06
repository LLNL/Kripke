#include<Kripke.h>
#include<Kripke/PartitionSpace.h>
#include<array>

using namespace Kripke;

PartitionSpace::PartitionSpace(Kripke::Comm &base_comm,
  size_t P, size_t Q, size_t Rx, size_t Ry, size_t Rz) :
  m_comm_all(base_comm),
  m_local_num_sdom{{0,0,0,0,0,0,0}},
  m_global_num_sdom{{0,0,0,0,0,0,0}},
  m_global_sdom_lower{{0,0,0,0,0,0,0}},
  m_proc_layout(P, Q, Rx, Ry, Rz),
  m_proc_xyz_layout(Rx, Ry, Rz),
  m_local_sdom_layout(0,0,0,0,0),
  m_local_sdom_xyz_layout(0,0,0),
  m_global_sdom_layout(0,0,0,0,0)
{
  size_t num_ranks = P*Q*Rx*Ry*Rz;

  // Check that our number of ranks is compatible
  KRIPKE_ASSERT(num_ranks == base_comm.size(),
    "Number of MPI ranks must match decomposition, expected %lu ranks\n",
    (unsigned long)num_ranks);

  // Assign communicators for P,Q,R
  m_comm_space[SPACE_PQR] = base_comm;
  m_comm_space[SPACE_R] = base_comm;

  // Compute our rank in pqxyz space
  std::array<int, 5> rank{{0,0,0,0,0}};
  m_proc_layout.toIndices(base_comm.rank(),
      rank[0], rank[1], rank[2], rank[3], rank[4]);

  // Project out dimensions to get our rank coloring in x,y,z
  for(size_t space = 0;space < 5;++ space){

    // Project out space
    std::array<int, 5> proj = rank;
    proj[space] = 0;

    // Get the coloring of this processor in space
    int color = m_proc_layout(proj[0], proj[1], proj[2], proj[3], proj[4]);

    // Split the communicator
    m_comm_space[space] = base_comm.split(color, rank[space]);
  };

  // Project out the R color and rank
  int rank_r = m_proc_layout(0, 0, rank[2], rank[3], rank[4]);
  int color_r = m_proc_layout(rank[0], rank[1], 0, 0, 0);
  
  // Split our R communicator
  m_comm_space[SPACE_R] = base_comm.split(color_r, rank_r);
}


void PartitionSpace::setup_createSubdomains(
    size_t SP, size_t SQ, size_t Sx, size_t Sy, size_t Sz){
  
  size_t num_sdom = SP * SQ * Sx * Sy * Sz;

  m_local_num_sdom[SPACE_PQR] = num_sdom;
  m_local_num_sdom[SPACE_P] = SP;
  m_local_num_sdom[SPACE_Q] = SQ;
  m_local_num_sdom[SPACE_RX] = Sx;
  m_local_num_sdom[SPACE_RY] = Sy;
  m_local_num_sdom[SPACE_RZ] = Sz;
  m_local_num_sdom[SPACE_R] = Sx * Sy * Sz;

  // Generate a local layout of subdomains
  m_local_sdom_layout = RAJA::Layout<5>(SP, SQ, Sx, Sy, Sz);

  // Create a local mapping from SR <==> (Sx, Sy, Sz)
  m_local_sdom_xyz_layout = RAJA::Layout<3>(Sx, Sy, Sz);

  // Compute global subdomain layout (and our local lower indices)
  for(size_t space = 0;space < NUM_SPACES;++ space){

    // Get the communicator for this space
    Kripke::Comm const &comm = m_comm_space[space];

    // Compute the total number of subdomains in this space's partition
    m_global_num_sdom[space] = comm.allReduceSumLong(m_local_num_sdom[space]);

    // Compute our lower offset into that global count
    m_global_sdom_lower[space] = comm.scanSumLong(m_local_num_sdom[space]) -
        m_local_num_sdom[space];

  }

  m_global_sdom_layout = RAJA::Layout<5>(m_global_num_sdom[SPACE_P],
                                         m_global_num_sdom[SPACE_Q],
                                         m_global_num_sdom[SPACE_RX],
                                         m_global_num_sdom[SPACE_RY],
                                         m_global_num_sdom[SPACE_RZ]);


}


size_t PartitionSpace::getNumSubdomains(Kripke::SPACE space) const{
  return m_local_num_sdom[space];
}

size_t PartitionSpace::getGlobalNumSubdomains(Kripke::SPACE space) const{
  return m_global_num_sdom[space];
}

size_t PartitionSpace::subdomainToSpace(
    Kripke::SPACE space, SdomId sdom_id) const
{
  // Map the subdomain id back to the bases spaces, P, Q, Rx, Ry, Rz
  std::array<int, 5> idx;
  m_local_sdom_layout.toIndices(*sdom_id,
      idx[0], idx[1], idx[2], idx[3], idx[4]);

  switch(space){
  // For the full space, it's just 1:1 with sdom's
  case SPACE_PQR:  return *sdom_id;

  // For R, we need to map back from Rx, Ry, Rz
  case SPACE_R:    return m_local_sdom_xyz_layout(idx[SPACE_RX],
                                                  idx[SPACE_RY],
                                                  idx[SPACE_RZ]);

  // For all other spaces, we just return their already computed index
  default:         return idx[space];
  }
}


SdomId PartitionSpace::spaceToSubdomain(
    Kripke::SPACE space, size_t sdom_space) const
{
  // build up indices in the P, Q, Rx, Ry, Rz space
  std::array<int, 5> idx{{0, 0, 0, 0, 0}};
  switch(space){
  // For the full space, it's 1:1 with sdom's
  case SPACE_PQR: return SdomId(sdom_space);

  // First, we need to compute Rx, Ry, and Rz
  case SPACE_R:   m_local_sdom_xyz_layout.toIndices(sdom_space,
                                                    idx[SPACE_RX],
                                                    idx[SPACE_RY],
                                                    idx[SPACE_RZ]);
                  break;

  // For basis spaces, just assign the index
  default:        idx[space] = sdom_space;
  }

  // convert those indices to a subdomain
  return SdomId{m_local_sdom_layout(idx[0], idx[1], idx[2], idx[3], idx[4])};
}


void PartitionSpace::print() const{
  if(m_comm_all.rank() == 0){
    printf("Decomposition Details:   Procs:      Subdomains (local/global):\n");
    printf("-----------------------  ----------  --------------------------\n");
    printf("  Total:                 %-10d  %d / %d\n",
        (int)m_comm_all.size(),
        (int)getNumSubdomains(),
        (int)(m_global_num_sdom[SPACE_P] *
              m_global_num_sdom[SPACE_Q] *
              m_global_num_sdom[SPACE_RX] *
              m_global_num_sdom[SPACE_RY] *
              m_global_num_sdom[SPACE_RZ]));
    printf("  (P) Energy:            %-10d  %d / %d\n",
        (int)m_comm_space[SPACE_P].size(),
        (int)m_local_num_sdom[SPACE_P],
        (int)m_global_num_sdom[SPACE_P]);
    printf("  (Q) Direction:         %-10d  %d / %d\n",
        (int)m_comm_space[SPACE_Q].size(),
        (int)m_local_num_sdom[SPACE_Q],
        (int)m_global_num_sdom[SPACE_Q]);
    printf("  (R) Space:             %-10d  %d / %d\n",
        (int)m_comm_space[SPACE_R].size(),
        (int)m_local_num_sdom[SPACE_R],
        (int)m_global_num_sdom[SPACE_R]);
    printf("  Space in XYZ:          %dx%dx%d       %dx%dx%d / %dx%dx%d\n",
        (int)m_comm_space[SPACE_RX].size(),
        (int)m_comm_space[SPACE_RY].size(),
        (int)m_comm_space[SPACE_RZ].size(),

        (int)m_local_num_sdom[SPACE_RX],
        (int)m_local_num_sdom[SPACE_RY],
        (int)m_local_num_sdom[SPACE_RZ],

        (int)m_global_num_sdom[SPACE_RX],
        (int)m_global_num_sdom[SPACE_RY],
        (int)m_global_num_sdom[SPACE_RZ]);

  }
}
