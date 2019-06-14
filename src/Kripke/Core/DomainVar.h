//
// Copyright (c) 2014-19, Lawrence Livermore National Security, LLC
// and Kripke project contributors. See the COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
//

#ifndef KRIPKE_CORE_DOMAIN_VAR_H__
#define KRIPKE_CORE_DOMAIN_VAR_H__

#include <Kripke.h>
#include <Kripke/Core/BaseVar.h>
#include <Kripke/Core/PartitionSpace.h>

namespace Kripke {
namespace Core {

  /**
   * Base class for variables that are defined over a PartitionSpace's
   * subdomains
   */
  class DomainVar : public Kripke::Core::BaseVar {
    public:
      DomainVar() = default;
      virtual ~DomainVar() = default;

      // Do not allow assignment or copy construction
      DomainVar(DomainVar const &) = delete;
      DomainVar& operator=(DomainVar const &) = delete;

      RAJA_INLINE
      size_t getNumSubdomains() const {
        return m_subdomain_to_chunk.size();
      }

      RAJA_INLINE
      std::vector<Kripke::SdomId> const &getWorkList() const {
        return m_work_list;
      }


      RAJA_INLINE
      void dump() const {
        printf("DomainVar:\n");

        printf("  m_subdomain_to_chunk: ");
        for(auto x : m_subdomain_to_chunk){printf("%lu ", (unsigned long)x);}
        printf("\n");

        printf("  m_chunk_to_subdomain: ");
        for(auto x : m_chunk_to_subdomain){printf("%lu ", (unsigned long)x);}
        printf("\n");

        printf("  m_work_list: ");
        for(auto x : m_work_list){printf("%d ", (int)*x);}
        printf("\n");
      }

    protected:

      void setup_initChunks(Kripke::Core::PartitionSpace const &pspace,
          Kripke::Core::SPACE space);

      void setup_initChunks(Kripke::Core::DomainVar const &clone_from);

      std::vector<size_t> m_subdomain_to_chunk;
      std::vector<size_t> m_chunk_to_subdomain;
      std::vector<Kripke::SdomId> m_work_list;
  };

} } // namespace

#endif
