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

#ifndef KRIPKE_CORE_PARTITION_SPACE_H__
#define KRIPKE_CORE_PARTITION_SPACE_H__

#include <Kripke.h>
#include <Kripke/Core/Comm.h>
#include <Kripke/Core/DataStore.h>
#include <array>

namespace Kripke {
namespace Core {

enum SPACE {
  SPACE_P = 0,
  SPACE_Q,
  SPACE_RX,
  SPACE_RY,
  SPACE_RZ,
  SPACE_R,
  SPACE_PR,
  SPACE_PQR,
  SPACE_NULL,
  NUM_SPACES
};

/**
 *  Defines a decomposition of the phase space by subdomains.
 */
class PartitionSpace : public Kripke::Core::BaseVar {
  public:
    using SdomCoord = std::array<ptrdiff_t, 5>;

    PartitionSpace(Kripke::Core::Comm &base_comm, 
      size_t P, size_t Q, size_t Rx, size_t Ry, size_t Rz);

    virtual ~PartitionSpace() = default;

    void setup_createSubdomains(
        size_t SP, size_t SQ, size_t Sx, size_t Sy, size_t Sz);

    void createSubdomainData(Kripke::Core::DataStore &data_store) const;

    size_t getNumSubdomains(Kripke::Core::SPACE space = SPACE_PQR) const;
    size_t getGlobalNumSubdomains(Kripke::Core::SPACE space = SPACE_PQR) const;

    SdomCoord sdomIdToCoord(Kripke::SdomId sdom_id) const;
    Kripke::SdomId coordToSdomId(SdomCoord coord) const;

    SdomCoord coordToGlobalCoord(SdomCoord local_coord) const;
    SdomCoord globalSdomIdToCoord(Kripke::GlobalSdomId global_sdom_id) const;
    Kripke::GlobalSdomId coordToGlobalSdomId(SdomCoord global_coord) const;



    Kripke::Core::Comm const &getComm(SPACE space) const {
      return m_comm_space[space];
    }
    
    size_t subdomainToSpace(Kripke::Core::SPACE space, SdomId sdom_id) const;
    SdomId spaceToSubdomain(Kripke::Core::SPACE space, size_t sdom_space) const;

    void print() const;

  private:
    Kripke::Core::Comm m_comm_all;

    // Parallel decomposition of comm_all
    Kripke::Core::Comm m_comm_space[NUM_SPACES];
    
    // Decomposition of ranks into subdomains
    std::array<long, NUM_SPACES> m_local_num_sdom;
    std::array<long, NUM_SPACES> m_global_num_sdom;
    std::array<long, NUM_SPACES> m_global_sdom_lower;


    std::array<RAJA::Layout<5>, NUM_SPACES> m_local_sdom_space_layout;

    RAJA::Layout<5> m_proc_layout;
    RAJA::Layout<3> m_proc_xyz_layout;
    RAJA::Layout<5> m_global_sdom_layout;

};



template<typename ELEMENT, typename ... IDX_TYPES>
class Field;

} // namespace Core

using Field_SdomId2GlobalSdomId = Kripke::Core::Field<GlobalSdomId, SdomId>;
using Field_GlobalSdomId2Rank = Kripke::Core::Field<long, GlobalSdomId>;
using Field_GlobalSdomId2SdomId = Kripke::Core::Field<SdomId, GlobalSdomId>;


} // namespace

#endif
