//
// Copyright (c) 2014-22, Lawrence Livermore National Security, LLC
// and Kripke project contributors. See the COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
//

#include <Kripke/Core/Set.h>


using namespace Kripke;
using namespace Kripke::Core;



/*****************************************************************************
 *
 *  Kripke::Core::Set
 *
 *****************************************************************************/


Set::Set() :
    m_global_size(0)
{

}


size_t Set::dimSize(Kripke::SdomId sdom_id, size_t ) const{
  return size(sdom_id);
}


/*****************************************************************************
 *
 *  Kripke::RangeSet
 *
 *****************************************************************************/

RangeSet::RangeSet(Kripke::Core::PartitionSpace const &pspace, Kripke::Core::SPACE space,
            std::vector<size_t> const &local_sizes) :
    m_space(space)
{
  setup_setupByLocalSize(pspace, local_sizes);
}


void RangeSet::setup_setupByLocalSize(Kripke::Core::PartitionSpace const &pspace,
    std::vector<size_t> const &local_sizes)
{

  Comm const &comm = pspace.getComm(m_space);

  // Figure out number of subdomains and chunks
  setup_initChunks(pspace, m_space);
  size_t num_chunks = m_chunk_to_subdomain.size();

  KRIPKE_ASSERT(local_sizes.size() == num_chunks,
    "Space %d has %lu subdomains, but provided %lu subdomains",
    (int)m_space,
    (unsigned long)m_chunk_to_subdomain.size(),
    (unsigned long)local_sizes.size());



  // Compute global size
  long total_local = 0;
  for(size_t s : local_sizes){
    total_local += s;
  }
  m_global_size = comm.allReduceSumLong(total_local);


  // Copy in local subdomain sizes
  m_chunk_to_size = local_sizes;


  // Compute global offsets for each chunk
  m_chunk_to_lower.resize(num_chunks);
  m_chunk_to_lower[0] = comm.scanSumLong(total_local) - total_local;
  for(size_t i = 1;i < num_chunks;++ i){
    m_chunk_to_lower[i] = m_chunk_to_lower[i-1] + m_chunk_to_size[i-1];
  }

}


/*****************************************************************************
 *
 *  Kripke::LocalRangeSet
 *
 *****************************************************************************/

LocalRangeSet::LocalRangeSet(Kripke::Core::PartitionSpace const &pspace,
            size_t local_size)
{

  Comm const &comm = pspace.getComm(SPACE_PQR);

  // Figure out number of subdomains and chunks
  setup_initChunks(pspace, SPACE_NULL);
  size_t num_chunks = m_chunk_to_subdomain.size();

  KRIPKE_ASSERT(num_chunks == 1, "Something's wrong, SPACE_NULL should have 1");

  // Compute global size
  m_global_size = comm.allReduceSumLong(local_size);

  // Copy in local subdomain size
  m_chunk_to_size = {local_size};

  // Compute global offsets for each chunk
  m_chunk_to_lower.resize(num_chunks);
  m_chunk_to_lower[0] = comm.scanSumLong(local_size) - local_size;
}

/*****************************************************************************
 *
 *  Kripke::GlobalRangeSet
 *
 *****************************************************************************/
GlobalRangeSet::GlobalRangeSet(Kripke::Core::PartitionSpace const &pspace,
    size_t global_size)
{
  setup_setGlobalSize(pspace, global_size);
}


GlobalRangeSet::GlobalRangeSet(Kripke::Core::PartitionSpace const &pspace, Kripke::Core::Set &parent_set)
{
  setup_setGlobalSize(pspace, parent_set.globalSize());
}


void GlobalRangeSet::setup_setGlobalSize(Kripke::Core::PartitionSpace const &pspace,
    size_t global_size)
{

  setup_initChunks(pspace, SPACE_NULL);

  // Kripke::Core::Set
  m_chunk_to_size.resize(1, global_size);
  m_chunk_to_lower.resize(1, 0);
  m_global_size = global_size;
}


