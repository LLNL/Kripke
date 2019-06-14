//
// Copyright (c) 2014-19, Lawrence Livermore National Security, LLC
// and Kripke project contributors. See the COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
//


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
