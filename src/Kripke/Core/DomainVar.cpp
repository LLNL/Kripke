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


#include <Kripke/Core/DomainVar.h>

using namespace Kripke;
using namespace Kripke::Core;


void DomainVar::setup_initChunks(Kripke::Core::PartitionSpace const &pspace,
          Kripke::Core::SPACE space)
{

  size_t num_subdomains = pspace.getNumSubdomains();
  size_t num_chunks = pspace.getNumSubdomains(space);


  // Map subdomains to chunks
  m_subdomain_to_chunk.resize(num_subdomains);
  for(SdomId sdom_id{0};sdom_id < (int)num_subdomains;++ sdom_id){
    size_t chunk_id = pspace.subdomainToSpace(space, sdom_id);
    m_subdomain_to_chunk[*sdom_id] = chunk_id;
  }

  // Map chunks to subdomains
  m_chunk_to_subdomain.resize(num_chunks);
  m_work_list.resize(num_chunks);
  for(size_t chunk_id = 0;chunk_id < num_chunks;++ chunk_id){
    SdomId sdom_id = pspace.spaceToSubdomain(space, chunk_id);
    m_chunk_to_subdomain[chunk_id] = *sdom_id;
    m_work_list[chunk_id] = sdom_id;
  }

}


void DomainVar::setup_initChunks(Kripke::Core::DomainVar const &clone_from){
  m_subdomain_to_chunk = clone_from.m_subdomain_to_chunk;
  m_chunk_to_subdomain = clone_from.m_chunk_to_subdomain;
  m_work_list = clone_from.m_work_list;
}
