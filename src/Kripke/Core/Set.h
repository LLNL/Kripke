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

#ifndef KRIPKE_CORE_SET_H__
#define KRIPKE_CORE_SET_H__

#include <Kripke.h>
#include <Kripke/Core/DomainVar.h>
#include <Kripke/Core/PartitionSpace.h>
#include <map>
#include <vector>

namespace Kripke {

namespace Core {
  /**
   * Base class for defining an ordered set used for dimensioning Fields
   */
  class Set : public Kripke::Core::DomainVar {
    public:
      Set();
      virtual ~Set() = default;

      // Don't allow copy construction
      Set(Set const &) = delete;

      /**
       * Returns the number of elements in this subdomain.
       */
      RAJA_INLINE
      size_t size(Kripke::SdomId sdom_id) const {
        size_t chunk_id = m_subdomain_to_chunk[*sdom_id];
        return m_chunk_to_size[chunk_id];
      }


      /**
       * Returns the range of a subdomain using a RAJA::RangeSegment
       */
      RAJA_INLINE
      RAJA::RangeSegment range(Kripke::SdomId sdom_id) const {
        return RAJA::RangeSegment(0, size(sdom_id));
      }


      /**
       * Returns the first global index for this subdomain.
       */
      RAJA_INLINE
      size_t lower(Kripke::SdomId sdom_id) const {
        size_t chunk_id = m_subdomain_to_chunk[*sdom_id];
        return m_chunk_to_lower[chunk_id];
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
      virtual size_t getNumDimensions() const = 0;

      /**
       * Returns the size of this Set along the specified dimension
       */
      virtual size_t dimSize(Kripke::SdomId sdom_id, size_t dim) const;

    protected:
      std::vector<size_t> m_chunk_to_size;
      std::vector<size_t> m_chunk_to_lower;
      size_t m_global_size;
  };


  class RangeSet : public Kripke::Core::Set {
    public:
      RangeSet(Kripke::Core::PartitionSpace const &pspace, Kripke::Core::SPACE space,
          std::vector<size_t> const &local_sizes);

      virtual ~RangeSet() = default;

      RAJA_INLINE
      virtual size_t getNumDimensions() const{return 1;}

    private:
      void setup_setupByLocalSize(Kripke::Core::PartitionSpace const &pspace,
          std::vector<size_t> const &local_sizes);
      Kripke::Core::SPACE m_space;
  };


  class LocalRangeSet : public Kripke::Core::Set {
    public:
      LocalRangeSet(Kripke::Core::PartitionSpace const &pspace, size_t local_size);

      virtual ~LocalRangeSet() = default;

      RAJA_INLINE
      virtual size_t getNumDimensions() const{return 1;}
  };


  class GlobalRangeSet : public Kripke::Core::Set {
    public:
      GlobalRangeSet(Kripke::Core::PartitionSpace const &pspace, size_t global_size);
      GlobalRangeSet(Kripke::Core::PartitionSpace const &pspace, Kripke::Core::Set &parent_set);

      virtual ~GlobalRangeSet() = default;

      RAJA_INLINE
      virtual size_t getNumDimensions() const{return 1;}

    private:
      void setup_setGlobalSize(Kripke::Core::PartitionSpace const &pspace, size_t global_size);
  };


  template<size_t NUM_SETS>
  class ProductSet : public Kripke::Core::Set {
    public:

      using LayoutType = RAJA::Layout<NUM_SETS>;

      template<typename ... SPAN>
      ProductSet(Kripke::Core::PartitionSpace &pspace, Kripke::Core::SPACE space,
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

      /**
       * Returns the size of this Set along the specified dimension
       */
      virtual size_t dimSize(Kripke::SdomId sdom_id, size_t dim) const{
        return m_spanned_sets[dim]->size(sdom_id);
      }

      RAJA_INLINE
      LayoutType getLayout(Kripke::SdomId sdom_id) const {

        std::array<RAJA::Index_type, NUM_SETS> sizes;
        for(size_t dim = 0;dim < NUM_SETS;++ dim){
          sizes[dim] = dimSize(sdom_id, dim);
        }

        //auto perm = camp::make_idx_seq<NUM_SETS>::array();
				auto perm = RAJA::as_array<RAJA::MakePerm<NUM_SETS>>::get();

        return RAJA::make_permuted_layout(sizes, perm);
      }

    private:

      /**
       * Helper function to expand variadic arguments to the constructor.
       */
      void setup_setSpannedSets(
          std::array<Kripke::Core::Set const *, NUM_SETS> const &spanned_sets){
        m_spanned_sets = spanned_sets;

        size_t num_chunks = m_chunk_to_subdomain.size();
        m_chunk_to_size.resize(num_chunks, 1);
        m_chunk_to_lower.resize(num_chunks, 0);
        for(size_t chunk_id = 0;chunk_id < num_chunks;++ chunk_id){
          Kripke::SdomId sdom_id(m_chunk_to_subdomain[chunk_id]);
          for(size_t set_id = 0;set_id < NUM_SETS;++ set_id){
            m_chunk_to_size[chunk_id] *= spanned_sets[set_id]->size(sdom_id);
          }
        }

        // Compute global size
        m_global_size = 1;
        for(size_t set_id = 0;set_id < NUM_SETS;++ set_id){
          m_global_size *= spanned_sets[set_id]->globalSize();
        }
      }

      static const size_t s_num_sets = NUM_SETS;

      std::array<Kripke::Core::Set const *, NUM_SETS> m_spanned_sets;


  };



} } // namespace

#endif

