//
// Copyright (c) 2014-19, Lawrence Livermore National Security, LLC
// and Kripke project contributors. See the COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
//

#ifndef KRIPKE_CORE_FIELD_H__
#define KRIPKE_CORE_FIELD_H__

#include <Kripke.h>
#include <Kripke/Core/VarLayout.h>
#include <Kripke/Core/DataStore.h>
#include <Kripke/Core/DomainVar.h>
#include <Kripke/Core/Set.h>
#include <vector>

#ifdef KRIPKE_USE_CHAI
#define DEBUG
#include <chai/ManagedArray.hpp>
#undef DEBUG
#endif

namespace Kripke {
namespace Core {
  /**
   * Base class for Field which provides storage allocation
   */
  template<typename ELEMENT>
  class FieldStorage : public Kripke::Core::DomainVar {
    public:
      using ElementType = ELEMENT;

#ifndef KRIPKE_USE_CHAI
      using ElementPtr = ELEMENT*;
#else
      using ElementPtr = chai::ManagedArray<ELEMENT>;
#endif

      using Layout1dType = RAJA::TypedLayout<RAJA::Index_type, camp::tuple<RAJA::Index_type>>;
      using View1dType = RAJA::View<ElementType, Layout1dType, ElementPtr>;


      explicit FieldStorage(Kripke::Core::Set const &spanned_set) :
        m_set(&spanned_set)
      {

        // initialize our decomposition to match that of the specified set
        setup_initChunks(spanned_set);

        // allocate all of our chunks, and create layouts for each one
        size_t num_chunks = m_chunk_to_subdomain.size();
        m_chunk_to_size.resize(num_chunks, 0);
#ifndef KRIPKE_USE_CHAI
        m_chunk_to_data.resize(num_chunks, nullptr);
#else
        m_chunk_to_data.resize(num_chunks);
#endif

        for(size_t chunk_id = 0;chunk_id < num_chunks;++ chunk_id){

          // Get the size of the subdomain from the set
          SdomId sdom_id(m_chunk_to_subdomain[chunk_id]);
          size_t sdom_size = spanned_set.size(sdom_id);

          m_chunk_to_size[chunk_id] = sdom_size;
#ifndef KRIPKE_USE_CHAI
          m_chunk_to_data[chunk_id] = new ElementType[sdom_size];
#else
          m_chunk_to_data[chunk_id].allocate(sdom_size, chai::CPU,
              [=](const chai::PointerRecord* record, chai::Action action, chai::ExecutionSpace space){
                /*printf("CHAI[%s, %d]: ", BaseVar::getName().c_str(), (int)chunk_id);
                switch(action){
                case chai::ACTION_ALLOC: printf("ALLOC "); break;
                case chai::ACTION_FREE: printf("FREE  "); break;
                case chai::ACTION_MOVE: printf("MOVE  "); break;
                default: printf("UNKNOWN ");
                }

                switch(space){
                case chai::CPU: printf("CPU "); break;
#ifdef KRIPKE_USE_CUDA
                case chai::GPU: printf("GPU  "); break;
#endif
                default: printf("UNK ");
                }

                printf("%lu bytes\n", (unsigned long) bytes);
*/
              }

          );
#endif
        }
      }

      virtual ~FieldStorage(){
#ifndef KRIPKE_USE_CHAI
        for(auto i : m_chunk_to_data){
          delete[] i;
        }
#endif
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
      ElementType *getData(Kripke::SdomId sdom_id) const {
        KRIPKE_ASSERT(*sdom_id < (int)m_subdomain_to_chunk.size(),
            "sdom_id(%d) >= num_subdomains(%d)",
            (int)*sdom_id,
            (int)(int)m_subdomain_to_chunk.size());
        size_t chunk_id = m_subdomain_to_chunk[*sdom_id];

#ifndef KRIPKE_USE_CHAI
        return  m_chunk_to_data[chunk_id];
#else
        // use pointer conversion to get host pointer
        ElementType *ptr = m_chunk_to_data[chunk_id];

        // return host pointer
        return(ptr);

#endif
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
#ifndef KRIPKE_USE_CHAI
      using ElementPtr = ELEMENT*;
#else
      using ElementPtr = chai::ManagedArray<ELEMENT>;
#endif

      static constexpr size_t NumDims = sizeof...(IDX_TYPES);

      using DefaultLayoutType = RAJA::TypedLayout<RAJA::Index_type, camp::tuple<IDX_TYPES...>>;
      using DefaultViewType = RAJA::View<ElementType, DefaultLayoutType, ElementPtr>;

      template<typename Order>
      Field(Kripke::Core::Set const &spanned_set, Order) :
        Parent(spanned_set)
      {

        KRIPKE_ASSERT(NumDims == spanned_set.getNumDimensions(),
            "Number of dimensions must match between Field<%d> and Set<%d>\n",
            (int)NumDims, (int)spanned_set.getNumDimensions());

        auto perm = LayoutInfo<Order, IDX_TYPES...>::getPermutation();

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


      template<typename Order>
      RAJA_INLINE
      auto getViewOrder(Kripke::SdomId sdom_id) const ->
        ViewType<Order, ElementType, ElementPtr, IDX_TYPES...>
      {
        size_t chunk_id = Parent::m_subdomain_to_chunk[*sdom_id];

        ElementPtr ptr = Parent::m_chunk_to_data[chunk_id];

        using LInfo = LayoutInfo<Order, IDX_TYPES...>;
        using LType = typename LInfo::Layout;

        LType layout = RAJA::make_stride_one<LInfo::stride_one_dim>(m_chunk_to_layout[chunk_id]);

        return ViewType<Order, ElementType, ElementPtr, IDX_TYPES...>(ptr, layout);
      }



      RAJA_INLINE
      void dump() const {
        printf("Field<>:\n");
        printf("  name:  %s\n", BaseVar::getName().c_str());
        printf("  m_set: %p\n", Parent::m_set);

        printf("  m_chunk_to_size: ");
        for(auto x : Parent::m_chunk_to_size){printf("%lu ", (unsigned long)x);}
        printf("\n");

#ifndef KRIPKE_USE_CHAI
        printf("  m_chunk_to_data: ");
        for(auto x : Parent::m_chunk_to_data){printf("%p ", x);}
        printf("\n");
#endif

        for(size_t chunk_id = 0;chunk_id < Parent::m_chunk_to_data.size();++ chunk_id){

          SdomId sdom_id(DomainVar::m_chunk_to_subdomain[chunk_id]);

          ElementType *ptr = Parent::getData(sdom_id);

          printf("Chunk %d Data: ", (int)chunk_id);
          for(size_t i = 0;i < Parent::m_chunk_to_size[chunk_id];++ i){
            printf(" %e", ptr[i]);
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

