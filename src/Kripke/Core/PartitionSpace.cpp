#include <Kripke/Core/PartitionSpace.h>

#include <Kripke.h>
#include <Kripke/Core/Field.h>
#include <Kripke/Kernel.h>
#include <Kripke/Core/Set.h>
#include <array>

using namespace Kripke;
using namespace Kripke::Core;

PartitionSpace::PartitionSpace(Kripke::Core::Comm &base_comm,
  size_t P, size_t Q, size_t Rx, size_t Ry, size_t Rz) :
  m_comm_all(base_comm),
  m_local_num_sdom{{0,0,0,0,0,0,0}},
  m_global_num_sdom{{0,0,0,0,0,0,0}},
  m_global_sdom_lower{{0,0,0,0,0,0,0}},
  m_proc_layout(P, Q, Rx, Ry, Rz),
  m_proc_xyz_layout(Rx, Ry, Rz),
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

  // Project out the PR color and rank
  int rank_pr = m_proc_layout(rank[0], 0, rank[2], rank[3], rank[4]);
  int color_pr = m_proc_layout(0, rank[1], 0, 0, 0);

  // Split our PR communicator
  m_comm_space[SPACE_PR] = base_comm.split(color_pr, rank_pr);
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
  m_local_num_sdom[SPACE_PR] = SP * Sx * Sy * Sz;
  m_local_num_sdom[SPACE_NULL] = 1;

  m_local_sdom_space_layout[SPACE_P] =    RAJA::Layout<5>(SP, 0, 0, 0, 0);
  m_local_sdom_space_layout[SPACE_Q] =    RAJA::Layout<5>(0, SQ, 0, 0, 0);
  m_local_sdom_space_layout[SPACE_RX] =   RAJA::Layout<5>(0, 0, Sx, 0, 0);
  m_local_sdom_space_layout[SPACE_RY] =   RAJA::Layout<5>(0, 0, 0, Sy, 0);
  m_local_sdom_space_layout[SPACE_RZ] =   RAJA::Layout<5>(0, 0, 0, 0, Sz);
  m_local_sdom_space_layout[SPACE_R] =    RAJA::Layout<5>(0, 0, Sx, Sy, Sz);
  m_local_sdom_space_layout[SPACE_PR] =   RAJA::Layout<5>(SP, 0, Sx, Sy, Sz);
  m_local_sdom_space_layout[SPACE_PQR] =  RAJA::Layout<5>(SP, SQ, Sx, Sy, Sz);
  m_local_sdom_space_layout[SPACE_NULL] = RAJA::Layout<5>(0, 0, 0, 0, 0);


  // Compute global subdomain layout (and our local lower indices)
  for(size_t space = 0;space < NUM_SPACES;++ space){

    // Get the communicator for this space
    Kripke::Core::Comm const &comm = m_comm_space[space];

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

/**
 * Creates Set and Field objects that describe the subdomain decomposition.
 * @param data_store The DataStore in which to create the objects
 */
void PartitionSpace::createSubdomainData(Kripke::Core::DataStore &data_store) const {

  PartitionSpace &pspace = data_store.getVariable<PartitionSpace>("pspace");

  // Create a Set that has exactly 1 element for each subdomain
  auto &set_sdomid = data_store.newVariable<LocalRangeSet>("Set/SdomId",
      *this,
      getNumSubdomains(SPACE_PQR));

  // Create a linearized version of the above
  auto &set_global_sdomid =
      data_store.newVariable<GlobalRangeSet>("Set/GlobalSdomIdLinear", pspace, set_sdomid);

  // Create a Field to store mappings from local subdomains to global
  auto &field_local_to_global =
      data_store.newVariable<Field_SdomId2GlobalSdomId>(
          "SdomId2GlobalSdomId", set_sdomid);

  auto &field_global_to_local =
       data_store.newVariable<Field_GlobalSdomId2SdomId>(
           "GlobalSdomId2SdomId", set_global_sdomid);
  Kripke::Kernel::kConst(field_global_to_local, SdomId{0});

  auto &field_global_to_rank =
       data_store.newVariable<Field_GlobalSdomId2Rank>(
           "GlobalSdomId2Rank", set_global_sdomid);
  Kripke::Kernel::kConst(field_global_to_rank, 0);


  size_t rank = m_comm_all.rank();

  for(SdomId sdom_id : set_sdomid.getWorkList()){
    auto local_to_global = field_local_to_global.getView(sdom_id);
    auto global_to_local = field_global_to_local.getView(sdom_id);
    auto global_to_rank = field_global_to_rank.getView(sdom_id);

    for(SdomId local{0};local < set_sdomid.size(sdom_id);++ local){
      //GlobalSdomId global(*local + set_sdomid.lower(sdom_id));

      // Get local subdomain coordinates
      SdomCoord local_coord = sdomIdToCoord(local);

      // Offset local coordinate to global coordinates
      SdomCoord global_coord = coordToGlobalCoord(local_coord);

      // Convert global coordinates to a GlobalSubdomainId
      GlobalSdomId global =  coordToGlobalSdomId(global_coord);

      local_to_global(local) = global;
      global_to_local(global) = local;
      global_to_rank(global) = rank;
    }

    // Perform collective to gather global addresses of all subdomains

    m_comm_all.allReduceSumLong(field_global_to_rank.getData(sdom_id),
                                field_global_to_rank.size(sdom_id));

    m_comm_all.allReduceSumInt((int*)field_global_to_local.getData(sdom_id),
                               field_global_to_local.size(sdom_id));
  }
}


size_t PartitionSpace::getNumSubdomains(Kripke::Core::SPACE space) const{
  return m_local_num_sdom[space];
}

size_t PartitionSpace::getGlobalNumSubdomains(Kripke::Core::SPACE space) const{
  return m_global_num_sdom[space];
}


PartitionSpace::SdomCoord PartitionSpace::sdomIdToCoord(Kripke::SdomId sdom_id) const{

  SdomCoord coord;

  m_local_sdom_space_layout[SPACE_PQR].toIndices(*sdom_id,
      coord[0], coord[1], coord[2], coord[3], coord[4]);

  return coord;
}
Kripke::SdomId PartitionSpace::coordToSdomId(SdomCoord coord) const{

  SdomId sdom_id(m_local_sdom_space_layout[SPACE_PQR](
      coord[0], coord[1], coord[2], coord[3], coord[4]));

  return sdom_id;
}

PartitionSpace::SdomCoord PartitionSpace::coordToGlobalCoord(SdomCoord local_coord) const{
  SdomCoord global_coord{{
    (ptrdiff_t)(local_coord[0] + m_global_sdom_lower[SPACE_P]),
    (ptrdiff_t)(local_coord[1] + m_global_sdom_lower[SPACE_Q]),
    (ptrdiff_t)(local_coord[2] + m_global_sdom_lower[SPACE_RX]),
    (ptrdiff_t)(local_coord[3] + m_global_sdom_lower[SPACE_RY]),
    (ptrdiff_t)(local_coord[4] + m_global_sdom_lower[SPACE_RZ]) }};

  return global_coord;
}
PartitionSpace::SdomCoord PartitionSpace::globalSdomIdToCoord(Kripke::GlobalSdomId global_sdom_id) const{

  SdomCoord coord;

  m_global_sdom_layout.toIndices(*global_sdom_id,
      coord[0], coord[1], coord[2], coord[3], coord[4]);

  return coord;
}
Kripke::GlobalSdomId PartitionSpace::coordToGlobalSdomId(SdomCoord global_coord) const{

  GlobalSdomId global_sdom_id(m_global_sdom_layout(
      global_coord[0], global_coord[1], global_coord[2], global_coord[3], global_coord[4]));

  return global_sdom_id;
}

size_t PartitionSpace::subdomainToSpace(
    Kripke::Core::SPACE space, SdomId sdom_id) const
{
  // Map the subdomain id back to the bases spaces, P, Q, Rx, Ry, Rz
  std::array<int, 5> idx;
  m_local_sdom_space_layout[SPACE_PQR].toIndices(*sdom_id,
      idx[0], idx[1], idx[2], idx[3], idx[4]);

  size_t space_id = m_local_sdom_space_layout[space](
      idx[0], idx[1], idx[2], idx[3], idx[4]);

  return space_id;
}


SdomId PartitionSpace::spaceToSubdomain(
    Kripke::Core::SPACE space, size_t space_id) const
{
  // build up indices in the P, Q, Rx, Ry, Rz space
  std::array<int, 5> idx{{0, 0, 0, 0, 0}};
  m_local_sdom_space_layout[space].toIndices(space_id,
      idx[0], idx[1], idx[2], idx[3], idx[4]);

  // convert those indices to a subdomain
  SdomId sdom_id{m_local_sdom_space_layout[SPACE_PQR](idx[0], idx[1], idx[2], idx[3], idx[4])};

  return sdom_id;
}


void PartitionSpace::print() const{
  if(m_comm_all.rank() == 0){
    printf("  Decomposition Space:   Procs:      Subdomains (local/global):\n");
    printf("  ---------------------  ----------  --------------------------\n");
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
    printf("  (Rx,Ry,Rz) R in XYZ:   %dx%dx%d       %dx%dx%d / %dx%dx%d\n",
        (int)m_comm_space[SPACE_RX].size(),
        (int)m_comm_space[SPACE_RY].size(),
        (int)m_comm_space[SPACE_RZ].size(),

        (int)m_local_num_sdom[SPACE_RX],
        (int)m_local_num_sdom[SPACE_RY],
        (int)m_local_num_sdom[SPACE_RZ],

        (int)m_global_num_sdom[SPACE_RX],
        (int)m_global_num_sdom[SPACE_RY],
        (int)m_global_num_sdom[SPACE_RZ]);

    printf("  (PQR) TOTAL:           %-10d  %d / %d\n",
        (int)m_comm_all.size(),
        (int)getNumSubdomains(),
        (int)(m_global_num_sdom[SPACE_P] *
              m_global_num_sdom[SPACE_Q] *
              m_global_num_sdom[SPACE_RX] *
              m_global_num_sdom[SPACE_RY] *
              m_global_num_sdom[SPACE_RZ]));

  }
}
