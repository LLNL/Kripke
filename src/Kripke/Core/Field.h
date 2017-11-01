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

#ifndef KRIPKE_CORE_FIELD_H__
#define KRIPKE_CORE_FIELD_H__

#include <Kripke.h>
#include <Kripke/Core/ArchLayout.h>
#include <Kripke/Core/DataStore.h>
#include <Kripke/Core/DomainVar.h>
#include <Kripke/Core/Set.h>
#include <vector>

namespace Kripke {
namespace Core {
  /**
   * Base class for Field which provides storage allocation
   */
  template<typename ELEMENT>
  class FieldStorage : public Kripke::Core::DomainVar {
    public:
      using ElementType = ELEMENT;
      using ElementPtr = ELEMENT*;

      using Layout1dType = Kripke::Core::LayoutType<AL_Default, RAJA::Index_type>;
      using View1dType = Kripke::Core::ViewType<AL_Default, ElementType, RAJA::Index_type>;


      explicit FieldStorage(Kripke::Core::Set const &spanned_set) :
        m_set(&spanned_set)
      {

        // initialize our decomposition to match that of the specified set
        setup_initChunks(spanned_set);

        // allocate all of our chunks, and create layouts for each one
        size_t num_chunks = m_chunk_to_subdomain.size();
        m_chunk_to_size.resize(num_chunks, 0);
        m_chunk_to_data.resize(num_chunks, nullptr);

        for(size_t chunk_id = 0;chunk_id < num_chunks;++ chunk_id){

          // Get the size of the subdomain from the set
          SdomId sdom_id(m_chunk_to_subdomain[chunk_id]);
          size_t sdom_size = spanned_set.size(sdom_id);

          m_chunk_to_size[chunk_id] = sdom_size;
          m_chunk_to_data[chunk_id] = new ElementType[sdom_size];
        }
      }

      virtual ~FieldStorage(){
        for(auto i : m_chunk_to_data){
          delete[] i;
        }
      }

      // Dissallow copy construction
      FieldStorage(FieldStorage<ElementType> const &) = delete;

      /**
       * Returns the number of elements in this subdomain.
       */
      RAJA_INLINE
      size_t size(Kripke::SdomId sdom_id) const {
        size_t chunk_id = m_subdomain_to_chunk[*sdom_id];
        return m_chunk_to_size[chunk_id];
      }


      RAJA_INLINE
      View1dType getView1d(Kripke::SdomId sdom_id) const {

        size_t chunk_id = m_subdomain_to_chunk[*sdom_id];

        ElementPtr ptr = m_chunk_to_data[chunk_id];
        size_t sdom_size = m_chunk_to_size[chunk_id];

        return View1dType(ptr, Layout1dType(sdom_size));
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
      Kripke::Core::Set const &getSet() const {
        return *m_set;
      }

    protected:
      Kripke::Core::Set const *m_set;
      std::vector<size_t> m_chunk_to_size;
      std::vector<ElementPtr> m_chunk_to_data;
  };

  /**
   * Defines a multi-dimensional data field defined over a Set
   */
  template<typename ELEMENT, typename ... IDX_TYPES>
  class Field : public Kripke::Core::FieldStorage<ELEMENT> {
    public:

      using Parent = Kripke::Core::FieldStorage<ELEMENT>;

      using ElementType = ELEMENT;
      using ElementPtr = ELEMENT*;

      static constexpr size_t NumDims = sizeof...(IDX_TYPES);

      using DefaultLayoutType = RAJA::TypedLayout<RAJA::Index_type, camp::tuple<IDX_TYPES...>>;
      using DefaultViewType = RAJA::View<ElementType, DefaultLayoutType>;


      Field(Kripke::Core::Set const &spanned_set) :
        Parent(spanned_set)
      {

        KRIPKE_ASSERT(NumDims == spanned_set.getNumDimensions(),
            "Number of dimensions must match between Field<%d> and Set<%d>\n",
            (int)NumDims, (int)spanned_set.getNumDimensions());

        // create layouts for each chunk
        size_t num_chunks = Parent::m_chunk_to_subdomain.size();
        m_chunk_to_layout.resize(num_chunks);
        for(size_t chunk_id = 0;chunk_id < num_chunks;++ chunk_id){

          // Create a layout using dim sizes from the Set, and permutation
          // defined by the layout function
          SdomId sdom_id(Parent::m_chunk_to_subdomain[chunk_id]);
          std::array<RAJA::Index_type, NumDims> sizes;
          for(size_t dim = 0;dim < NumDims;++ dim){
            sizes[dim] = spanned_set.dimSize(sdom_id, dim);
          }

          using LInfo = LayoutInfo<Layout_Default, IDX_TYPES...>;
          auto perm = LInfo::getPermutation();

          RAJA::Layout<NumDims, RAJA::Index_type> &layout =
              m_chunk_to_layout[chunk_id];
          layout = RAJA::make_permuted_layout<NumDims,RAJA::Index_type>(sizes, perm);
        }
      }

      virtual ~Field(){

      }



      RAJA_INLINE
      DefaultViewType getView(Kripke::SdomId sdom_id) const {

        size_t chunk_id = Parent::m_subdomain_to_chunk[*sdom_id];

        auto ptr = Parent::m_chunk_to_data[chunk_id];
        auto layout = m_chunk_to_layout[chunk_id];

        return DefaultViewType(ptr, layout);
      }


      template<typename AL>
      RAJA_INLINE
      auto getViewAL(Kripke::SdomId sdom_id) const ->
        ViewType<typename AL::LayoutFamily, ElementType, IDX_TYPES...>
      {
        size_t chunk_id = Parent::m_subdomain_to_chunk[*sdom_id];

        ElementPtr ptr = Parent::m_chunk_to_data[chunk_id];

        using LInfo = LayoutInfo<typename AL::LayoutFamily, IDX_TYPES...>;
        using LType = typename LInfo::Layout;

        LType layout = RAJA::make_stride_one<LInfo::stride_one_dim>(m_chunk_to_layout[chunk_id]);

        return ViewType<typename AL::LayoutFamily, ElementType, IDX_TYPES...>(ptr, layout);
      }



      RAJA_INLINE
      void dump() const {
        printf("Field<>:\n");
        printf("  m_set: %p\n", Parent::m_set);

        printf("  m_chunk_to_size: ");
        for(auto x : Parent::m_chunk_to_size){printf("%lu ", (unsigned long)x);}
        printf("\n");

        printf("  m_chunk_to_data: ");
        for(auto x : Parent::m_chunk_to_data){printf("%p ", x);}
        printf("\n");

        for(size_t chunk_id = 0;chunk_id < Parent::m_chunk_to_data.size();++ chunk_id){
          printf("Chunk %d Data: ", (int)chunk_id);
          for(size_t i = 0;i < Parent::m_chunk_to_size[chunk_id];++ i){
            //printf(" %d", RAJA::convertIndex<int>(Parent::m_chunk_to_data[chunk_id][i]));
            printf(" %e", Parent::m_chunk_to_data[chunk_id][i]);
          }
          printf("\n");
        }

        Kripke::Core::DomainVar::dump();
      }

    protected:
      std::vector<DefaultLayoutType> m_chunk_to_layout;
  };

} } // namespace

#endif

