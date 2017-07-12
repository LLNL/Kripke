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

#ifndef KRIPKE_FIELD_H__
#define KRIPKE_FIELD_H__

#include <Kripke.h>
#include <Kripke/DataStore.h>
#include <Kripke/DomainVar.h>
#include <Kripke/Set.h>
#include <vector>

namespace Kripke {

  /**
   * Base class for defining an ordered set used for dimensioning Fields
   */
  template<typename ELEMENT, typename ... IDX_TYPES>
  class Field : public Kripke::DomainVar {
    public:

      using ElementType = ELEMENT;
      using ElementPtr = ELEMENT*;

      static constexpr size_t NumDims = sizeof...(IDX_TYPES);

      using LayoutType = RAJA::TypedLayout<RAJA::Index_type, IDX_TYPES...>;
      //using LayoutType = RAJA::Layout<NumDims>;
      using ViewType = RAJA::View<ElementType, LayoutType>;

      using LayoutFunction = std::function<
                              std::array<size_t, NumDims>(Kripke::SdomId)
                             >;

      static std::array<size_t, NumDims> defaultLayout(Kripke::SdomId) {
        return VarOps::make_index_sequence<NumDims>::value;
      }

      Field(Kripke::Set const &spanned_set,
            LayoutFunction layout_fcn = defaultLayout) :
        m_set(&spanned_set)
      {

        KRIPKE_ASSERT(NumDims == spanned_set.getNumDimensions(),
            "Number of dimensions must match between Field and Set");

        // initialize our decomposition to match that of the specified set
        setup_initChunks(spanned_set);

        // allocate all of our chunks, and create layouts for each one
        size_t num_chunks = m_chunk_to_subdomain.size();
        m_chunk_to_size.resize(num_chunks, 0);
        m_chunk_to_data.resize(num_chunks, nullptr);
        m_chunk_to_layout.resize(num_chunks);
        for(size_t chunk_id = 0;chunk_id < num_chunks;++ chunk_id){

          // Get the size of the subdomain from the set
          SdomId sdom_id(m_chunk_to_subdomain[chunk_id]);
          size_t sdom_size = spanned_set.size(sdom_id);

          m_chunk_to_size[chunk_id] = sdom_size;
          m_chunk_to_data[chunk_id] = new ElementType[sdom_size];

          // Create a layout using dim sizes from the Set, and permutation
          // defined by the layout function
          std::array<int, NumDims> sizes;
          for(size_t dim = 0;dim < NumDims;++ dim){
            sizes[dim] = spanned_set.dimSize(sdom_id, dim);
          }
          auto perm = layout_fcn(sdom_id);

          RAJA::Layout<NumDims, RAJA::Index_type> &layout =
              m_chunk_to_layout[chunk_id];
          layout = RAJA::make_permuted_layout<NumDims,RAJA::Index_type>(sizes, perm);
        }
      }

      virtual ~Field(){
        for(auto i : m_chunk_to_data){
          delete[] i;
        }
      }


      /**
       * Returns the number of elements in this subdomain.
       */
      RAJA_INLINE
      size_t size(Kripke::SdomId sdom_id) const {
        size_t chunk_id = m_subdomain_to_chunk[*sdom_id];
        return m_chunk_to_size[chunk_id];
      }

      RAJA_INLINE
      ViewType getView(Kripke::SdomId sdom_id) const {

        size_t chunk_id = m_subdomain_to_chunk[*sdom_id];

        ElementPtr ptr = m_chunk_to_data[chunk_id];
        LayoutType layout = m_chunk_to_layout[chunk_id];

        return ViewType(ptr, layout);
      }

      RAJA_INLINE
      ElementPtr getData(Kripke::SdomId sdom_id) const {

        KRIPKE_ASSERT(*sdom_id < (int)m_subdomain_to_chunk.size(),
            "sdom_id(%d) >= num_subdomains(%d)",
            (int)*sdom_id,
            (int)(int)m_subdomain_to_chunk.size());
        size_t chunk_id = m_subdomain_to_chunk[*sdom_id];

        return  m_chunk_to_data[chunk_id];
      }


      RAJA_INLINE
      void dump() const {
        printf("Field<>:\n");
        printf("  m_set: %p\n", m_set);

        printf("  m_chunk_to_size: ");
        for(auto x : m_chunk_to_size){printf("%lu ", (unsigned long)x);}
        printf("\n");

        printf("  m_chunk_to_data: ");
        for(auto x : m_chunk_to_data){printf("%p ", x);}
        printf("\n");

        DomainVar::dump();
      }

    protected:
      Kripke::Set const *m_set;
      std::vector<size_t> m_chunk_to_size;
      std::vector<ElementPtr> m_chunk_to_data;
      std::vector<LayoutType> m_chunk_to_layout;
  };

} // namespace

#endif

