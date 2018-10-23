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
