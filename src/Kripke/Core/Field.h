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

#if defined(KRIPKE_USE_ZFP)
#include "zfparray1.h"
#endif

namespace Kripke {
namespace Core {

#if defined(KRIPKE_USE_ZFP)

  template <typename CONFIG, typename ELEMENT>
  struct BasicStorageTypeHelper {
    using config_type = CONFIG;
    using element_type = ELEMENT;

    using storage_type = void; // only used with ZFP
    using storage_pointer = element_type *;
    using element_pointer = element_type *;
    using element_reference = element_type &;

    static inline storage_pointer alloc(size_t size) {
      return new element_type[size];
    }

    static inline void free(storage_pointer store) {
      delete [] store;
    }

    static inline element_pointer get_data(storage_pointer store) {
      return store;
    }
  };

  struct field_storage_config {};

  template<typename T, typename = void>
  struct has_type : std::false_type { };

  template<typename T>
  struct has_type<T, decltype(sizeof(typename T::type), void())> : std::true_type { };

  template<typename T, typename = void>
  struct has_zfp_rate : std::false_type { };

  template<typename T>
  struct has_zfp_rate<T, decltype(std::declval<T>().zfp_rate, void())> : std::true_type { };

  template <
    typename ELEMENT,
    unsigned int N,
    bool = std::conditional< std::is_base_of<field_storage_config, ELEMENT>::value, has_type<ELEMENT>, std::false_type >::type::value,
    bool = has_zfp_rate<ELEMENT>::value
  >
  struct StorageTypeHelper;

  template <typename ELEMENT, unsigned int N>
  struct StorageTypeHelper<ELEMENT, N, false, false> : BasicStorageTypeHelper<ELEMENT, ELEMENT> {};

  template <typename ELEMENT, unsigned int N>
  struct StorageTypeHelper<ELEMENT, N, true, false> : BasicStorageTypeHelper<ELEMENT, typename ELEMENT::type> {};

  template <typename ELEMENT, unsigned int N>
  struct StorageTypeHelper<ELEMENT, N, true, true> {
    using config_type = ELEMENT;

    using element_type = typename config_type::type;

    using storage_type = zfp::array1<element_type>;
    using storage_pointer = storage_type *;
    using element_pointer = typename storage_type::pointer;
    using element_reference = typename storage_type::reference;

    static inline storage_pointer alloc(size_t size) {
      return new storage_type(size, config_type::zfp_rate);
    }

    static inline void free(storage_pointer store) {
      delete store;
    }

    static inline element_pointer get_data(storage_pointer store) {
      return element_pointer(store, 0);
    }
  };

#endif

  /**
   * Base class for Field which provides storage allocation
   */
#if defined(KRIPKE_USE_ZFP)
  template<typename ELEMENT, unsigned int N>
#else
  template<typename ELEMENT>
#endif
  class FieldStorage : public Kripke::Core::DomainVar {
    public:

#if defined(KRIPKE_USE_ZFP)
      using StorageHelper = StorageTypeHelper<ELEMENT, N>;
      using ConfigType = typename StorageHelper::config_type;
      using StoragePtr = typename StorageHelper::storage_pointer;
      using ElementType = typename StorageHelper::element_type;
      using ElementPtr = typename StorageHelper::element_pointer;
      using ElementRef = typename StorageHelper::element_reference;
#else
      using ElementType = ELEMENT;
      using StoragePtr = ElementType *;
      using ElementRef = ElementType;
#if defined(KRIPKE_USE_CHAI)
      using ElementPtr = chai::ManagedArray<ElementType>;
#else
      using ElementPtr = ElementType *;
#endif
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
#if defined(KRIPKE_USE_CHAI)
        m_chunk_to_data.resize(num_chunks);
#else
        m_chunk_to_data.resize(num_chunks, nullptr);
#endif

        for(size_t chunk_id = 0;chunk_id < num_chunks;++ chunk_id){

          // Get the size of the subdomain from the set
          SdomId sdom_id(m_chunk_to_subdomain[chunk_id]);
          size_t sdom_size = spanned_set.size(sdom_id);

          m_chunk_to_size[chunk_id] = sdom_size;
#if defined(KRIPKE_USE_CHAI)
          m_chunk_to_data[chunk_id].allocate(sdom_size, chai::CPU,
              [=](chai::Action action, chai::ExecutionSpace space, size_t bytes){
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
#elif defined(KRIPKE_USE_ZFP)
          m_chunk_to_data[chunk_id] = StorageHelper::alloc(sdom_size);
#else
          m_chunk_to_data[chunk_id] = new ElementType[sdom_size];
#endif
        }
      }

      virtual ~FieldStorage(){
#if defined(KRIPKE_USE_CHAI)
#elif defined(KRIPKE_USE_ZFP)
        for(auto i : m_chunk_to_data){
          StorageHelper::free(i);
        }
#else
        for(auto i : m_chunk_to_data){
          delete[] i;
        }
#endif
      }

      // Dissallow copy construction
#if defined(KRIPKE_USE_ZFP)
      FieldStorage(FieldStorage<ConfigType, N> const &) = delete;
#else
      FieldStorage(FieldStorage<ElementType> const &) = delete;
#endif

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

#if defined(KRIPKE_USE_ZFP)
        ElementPtr ptr = StorageHelper::get_data(m_chunk_to_data[chunk_id]);
#else
        ElementPtr ptr = m_chunk_to_data[chunk_id];
#endif
        size_t sdom_size = m_chunk_to_size[chunk_id];

        return View1dType(ptr, Layout1dType(sdom_size));
      }

      RAJA_INLINE
      StoragePtr getData(Kripke::SdomId sdom_id) const {
        KRIPKE_ASSERT(*sdom_id < (int)m_subdomain_to_chunk.size(),
            "sdom_id(%d) >= num_subdomains(%d)",
            (int)*sdom_id,
            (int)(int)m_subdomain_to_chunk.size());
        size_t chunk_id = m_subdomain_to_chunk[*sdom_id];


#if defined(KRIPKE_USE_CHAI)
        // use pointer conversion to get host pointer
        ElementType *ptr = m_chunk_to_data[chunk_id];
        // return host pointer
        return(ptr);
#else
        return  m_chunk_to_data[chunk_id];
#endif
      }


      RAJA_INLINE
      Kripke::Core::Set const &getSet() const {
        return *m_set;
      }

    protected:
      Kripke::Core::Set const *m_set;
      std::vector<size_t> m_chunk_to_size;
#if defined(KRIPKE_USE_ZFP)
      std::vector<StoragePtr> m_chunk_to_data;
#else
      std::vector<ElementPtr> m_chunk_to_data;
#endif
  };

  /**
   * Defines a multi-dimensional data field defined over a Set
   */
  template<typename ELEMENT, typename ... IDX_TYPES>
#if defined(KRIPKE_USE_ZFP)
  class Field : public Kripke::Core::FieldStorage<ELEMENT, sizeof...(IDX_TYPES)> {
#else
  class Field : public Kripke::Core::FieldStorage<ELEMENT> {
#endif
    public:

#if defined(KRIPKE_USE_ZFP)
      using Parent = Kripke::Core::FieldStorage<ELEMENT, sizeof...(IDX_TYPES)>;
#else
      using Parent = Kripke::Core::FieldStorage<ELEMENT>;
#endif

#if defined(KRIPKE_USE_ZFP)
      using StorageHelper = typename Parent::StorageHelper;
      using ElementType = typename Parent::ElementType;
      using ConfigType = typename Parent::ConfigType;
      using StoragePtr = typename Parent::StoragePtr;
      using ElementPtr = typename Parent::ElementPtr;
      using ElementRef = typename Parent::ElementRef;
#else
      using ElementType = typename Parent::ElementType;
      using ElementPtr = typename Parent::ElementPtr;
      using ElementRef = typename Parent::ElementRef;
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

#if defined(KRIPKE_USE_ZFP)
        ElementPtr ptr = StorageHelper::get_data(Parent::m_chunk_to_data[chunk_id]);
#else
        ElementPtr ptr = Parent::m_chunk_to_data[chunk_id];
#endif
        auto layout = m_chunk_to_layout[chunk_id];

        return DefaultViewType(ptr, layout);
      }


      template<typename Order>
      RAJA_INLINE
      auto getViewOrder(Kripke::SdomId sdom_id) const ->
        ViewType<Order, ElementType, ElementPtr, IDX_TYPES...>
      {
        size_t chunk_id = Parent::m_subdomain_to_chunk[*sdom_id];

#if defined(KRIPKE_USE_ZFP)
        ElementPtr ptr = StorageHelper::get_data(Parent::m_chunk_to_data[chunk_id]);
#else
        ElementPtr ptr = Parent::m_chunk_to_data[chunk_id];
#endif

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

#if !defined(KRIPKE_USE_CHAI)
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

