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

#ifndef KRIPKE_SET_H__
#define KRIPKE_SET_H__

#include <Kripke.h>
#include <Kripke/DomainVar.h>
#include <Kripke/PartitionSpace.h>
#include <map>
#include <vector>

namespace Kripke {

  /**
   * Base class for defining an ordered set used for dimensioning Fields
   */
  class Set : public Kripke::DomainVar {
    public:
      Set();
      virtual ~Set() = default;

      /**
       * Returns the number of elements in this subdomain.
       */
      RAJA_INLINE
      size_t size(Kripke::SdomId sdom_id) const {
        size_t chunk_id = m_subdomain_to_chunk[*sdom_id];
        return m_chunk_to_size[chunk_id];
      }

      /**
       * Returns the global number of unique elements in this set.
       */
      RAJA_INLINE
      size_t globalSize() const {
        return m_global_size;
      }

      /**
       * Returns the dimensionality of this Set.
       */
      virtual size_t getNumDimensions() const;

      /**
       * Returns the size of this Set along the specified dimension
       */
      virtual size_t dimSize(size_t dim, Kripke::SdomId sdom_id) const;

    protected:
      std::vector<size_t> m_chunk_to_size;
      std::vector<size_t> m_chunk_to_lower;
      size_t m_global_size;
  };


  class RangeSet : public Kripke::Set {
    public:
      RangeSet(Kripke::PartitionSpace &pspace, Kripke::SPACE space,
          std::vector<size_t> const &local_sizes);

      virtual ~RangeSet() = default;

    private:
      void setup_setupByLocalSize(std::vector<size_t> const &local_sizes);
      Kripke::PartitionSpace *m_pspace;
      Kripke::SPACE m_space;
  };


  class GlobalRangeSet : public Kripke::Set {
    public:
      GlobalRangeSet(Kripke::PartitionSpace &pspace, size_t global_size);
      explicit GlobalRangeSet(Kripke::Set &parent_set);

      virtual ~GlobalRangeSet() = default;

    private:
      void setup_setGlobalSize(size_t num_subdomains, size_t global_size);
  };


  template<size_t NUM_SETS>
  class ProductSet : public Kripke::Set {
    public:

      template<typename ... SPAN>
      ProductSet(Kripke::PartitionSpace &pspace, Kripke::SPACE space,
          SPAN const &... spanned_sets){
        static_assert(sizeof...(SPAN) == NUM_SETS,
            "Must provide same number of sets as dimensionality of ProductSet");

        setup_initChunks(pspace, space);
        setup_setSpannedSets({{(&spanned_sets)...}});

      }

      virtual ~ProductSet() = default;

      virtual size_t getNumDimensions() const{
        return s_num_sets;
      }

    private:

      /**
       * Helper function to expand variadic arguments to the constructor.
       */
      void setup_setSpannedSets(
          std::array<Kripke::BaseVar const *, NUM_SETS> const &spanned_sets){
        m_spanned_sets = spanned_sets;
      }

      static const size_t s_num_sets = NUM_SETS;

      std::array<Kripke::BaseVar const *, NUM_SETS> m_spanned_sets;


  };



} // namespace

#endif

