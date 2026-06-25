//
// Copyright (c) 2014-25, Lawrence Livermore National Security, LLC
// and Kripke project contributors. See the Kripke/COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
//

#ifndef KRIPKE_CORE_FIELD_H__
#define KRIPKE_CORE_FIELD_H__

#include <Kripke.h>
#include <Kripke/Core/VarLayout.h>
#include <Kripke/Core/DataStore.h>
#include <Kripke/Core/DomainVar.h>
#include <Kripke/Core/MemoryManager.h>
#include <Kripke/Core/Set.h>
#include <cstring>
#include <vector>

#ifdef KRIPKE_USE_CHAI
#define DEBUG
#include <umpire/Umpire.hpp>
#undef DEBUG
#endif

namespace Kripke {
namespace Core {
  template<typename ELEMENT, bool HOST_RESIDENT_NORMAL_CHAI_GPU, typename ... IDX_TYPES>
  class FieldWithPolicy;

  /**
   * Base class for Field which provides storage allocation
   */
  template<typename ELEMENT>
  class FieldStorage : public Kripke::Core::DomainVar {
    public:
      using ElementType = ELEMENT;
      using ElementPtr = ELEMENT*;

      using Layout1dType = RAJA::TypedLayout<RAJA::Index_type, camp::tuple<RAJA::Index_type>>;
      //FIXME: Remove the internal namespace when new RAJA release is out
      using View1dType = RAJA::internal::ViewBase<ElementType, ElementPtr, Layout1dType>;


      explicit FieldStorage(Kripke::Core::Set const &spanned_set
#ifdef KRIPKE_USE_CHAI
          , Kripke::ExecutionSpace allocation_space = Kripke::CPU
#endif
          ) :
        m_set(&spanned_set)
#ifdef KRIPKE_USE_CHAI
        , m_allocation_space(allocation_space)
#endif
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
          size_t bytes = sdom_size*sizeof(ElementType);

          m_chunk_to_size[chunk_id] = sdom_size;
#ifndef KRIPKE_USE_CHAI
          m_chunk_to_data[chunk_id] = new ElementType[sdom_size];
#else
          auto &chunk = m_chunk_to_data[chunk_id];
          if(m_allocation_space == Kripke::GPU){
            chunk.device_ptr = allocateDeviceBuffer(bytes);
            chunk.last_space = Kripke::GPU;
          }
          else{
            chunk.host_ptr = allocateHostBuffer(bytes);
            chunk.last_space = Kripke::CPU;
          }
#endif
        }
      }

      virtual ~FieldStorage(){
#ifndef KRIPKE_USE_CHAI
        for(auto i : m_chunk_to_data){
          delete[] i;
        }
#else
        auto host_allocator = MemoryManager::getHostAllocator();
        for(auto &chunk : m_chunk_to_data){
          if(chunk.device_ptr != nullptr){
            MemoryManager::getDeviceAllocator().deallocate(chunk.device_ptr);
          }
          if(chunk.host_ptr != nullptr){
            host_allocator.deallocate(chunk.host_ptr);
          }
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
        ElementPtr ptr = getHostData(sdom_id);
        size_t chunk_id = m_subdomain_to_chunk[*sdom_id];
        size_t sdom_size = m_chunk_to_size[chunk_id];

        return View1dType(ptr, Layout1dType(sdom_size));
      }

      RAJA_INLINE
      ElementType *getHostData(Kripke::SdomId sdom_id) const {
        KRIPKE_ASSERT(*sdom_id < (int)m_subdomain_to_chunk.size(),
            "sdom_id(%d) >= num_subdomains(%d)",
            (int)*sdom_id,
            (int)(int)m_subdomain_to_chunk.size());
        size_t chunk_id = m_subdomain_to_chunk[*sdom_id];

#ifndef KRIPKE_USE_CHAI
        return  m_chunk_to_data[chunk_id];
#else
        ensureHostCurrent(chunk_id);
        m_chunk_to_data[chunk_id].last_space = Kripke::CPU;
        return m_chunk_to_data[chunk_id].host_ptr;
#endif
      }

      RAJA_INLINE
      ElementType const *getHostDataConst(Kripke::SdomId sdom_id) const {
        KRIPKE_ASSERT(*sdom_id < (int)m_subdomain_to_chunk.size(),
            "sdom_id(%d) >= num_subdomains(%d)",
            (int)*sdom_id,
            (int)(int)m_subdomain_to_chunk.size());
        size_t chunk_id = m_subdomain_to_chunk[*sdom_id];

#ifndef KRIPKE_USE_CHAI
        return  m_chunk_to_data[chunk_id];
#else
        ensureHostCurrent(chunk_id);
        return m_chunk_to_data[chunk_id].host_ptr;
#endif
      }

      RAJA_INLINE
      ElementType *getDeviceData(Kripke::SdomId sdom_id) const {
        KRIPKE_ASSERT(*sdom_id < (int)m_subdomain_to_chunk.size(),
            "sdom_id(%d) >= num_subdomains(%d)",
            (int)*sdom_id,
            (int)(int)m_subdomain_to_chunk.size());
        size_t chunk_id = m_subdomain_to_chunk[*sdom_id];

#ifndef KRIPKE_USE_CHAI
        return  m_chunk_to_data[chunk_id];
#else
#if defined(KRIPKE_USE_CUDA) || defined(KRIPKE_USE_HIP)
        if(m_allocation_space == Kripke::GPU){
          ensureDeviceCurrent(chunk_id);
          m_chunk_to_data[chunk_id].last_space = Kripke::GPU;
          return m_chunk_to_data[chunk_id].device_ptr;
        }
        m_chunk_to_data[chunk_id].last_space = Kripke::CPU;
        return m_chunk_to_data[chunk_id].host_ptr;
#else
        ensureHostCurrent(chunk_id);
        m_chunk_to_data[chunk_id].last_space = Kripke::CPU;
        return m_chunk_to_data[chunk_id].host_ptr;
#endif
#endif
      }

      RAJA_INLINE
      ElementType *getData(Kripke::SdomId sdom_id) const {
        return getHostData(sdom_id);
      }

#ifdef KRIPKE_USE_CHAI
      RAJA_INLINE
      Kripke::ExecutionSpace getAllocationSpace() const {
        return m_allocation_space;
      }

      RAJA_INLINE
      void registerDeviceTouch(Kripke::SdomId sdom_id) {
#if defined(KRIPKE_USE_CUDA) || defined(KRIPKE_USE_HIP)
        if(m_allocation_space == Kripke::GPU){
          KRIPKE_ASSERT(*sdom_id < (int)m_subdomain_to_chunk.size(),
              "sdom_id(%d) >= num_subdomains(%d)",
              (int)*sdom_id,
              (int)(int)m_subdomain_to_chunk.size());
          size_t chunk_id = m_subdomain_to_chunk[*sdom_id];
          m_chunk_to_data[chunk_id].last_space = Kripke::GPU;
        }
#else
        (void)sdom_id;
#endif
      }
#endif


      RAJA_INLINE
      Kripke::Core::Set const &getSet() const {
        return *m_set;
      }

    protected:
#ifdef KRIPKE_USE_CHAI
      struct ChunkData {
        ElementType *host_ptr = nullptr;
        ElementType *device_ptr = nullptr;
        mutable Kripke::ExecutionSpace last_space = Kripke::CPU;
      };

      RAJA_INLINE
      static ElementType *allocateHostBuffer(size_t bytes) {
        if(bytes == 0){
          return nullptr;
        }
        return static_cast<ElementType*>(MemoryManager::getHostAllocator().allocate(bytes));
      }

      RAJA_INLINE
      static ElementType *allocateDeviceBuffer(size_t bytes) {
        if(bytes == 0){
          return nullptr;
        }
        return static_cast<ElementType*>(MemoryManager::getDeviceAllocator().allocate(bytes));
      }

      RAJA_INLINE
      void ensureHostCurrent(size_t chunk_id) const {
        if(m_allocation_space != Kripke::GPU ||
           m_chunk_to_data[chunk_id].device_ptr == nullptr ||
           m_chunk_to_data[chunk_id].last_space == Kripke::CPU){
          return;
        }

        if(m_chunk_to_data[chunk_id].host_ptr == nullptr){
          m_chunk_to_data[chunk_id].host_ptr = allocateHostBuffer(
              m_chunk_to_size[chunk_id]*sizeof(ElementType));
        }
        synchronizeDevice();
        MemoryManager::copy(m_chunk_to_data[chunk_id].host_ptr,
                            m_chunk_to_data[chunk_id].device_ptr,
                            m_chunk_to_size[chunk_id]*sizeof(ElementType));
        m_chunk_to_data[chunk_id].last_space = Kripke::CPU;
      }

      RAJA_INLINE
      void ensureDeviceCurrent(size_t chunk_id) const {
        if(m_allocation_space != Kripke::GPU ||
           m_chunk_to_data[chunk_id].last_space == Kripke::GPU){
          return;
        }

        if(m_chunk_to_data[chunk_id].device_ptr == nullptr){
          m_chunk_to_data[chunk_id].device_ptr = allocateDeviceBuffer(
              m_chunk_to_size[chunk_id]*sizeof(ElementType));
        }
        MemoryManager::copy(m_chunk_to_data[chunk_id].device_ptr,
                            m_chunk_to_data[chunk_id].host_ptr,
                            m_chunk_to_size[chunk_id]*sizeof(ElementType));
        m_chunk_to_data[chunk_id].last_space = Kripke::GPU;
      }

      RAJA_INLINE
      void synchronizeDevice() const {
#if defined(KRIPKE_USE_CUDA)
        RAJA::synchronize<RAJA::cuda_synchronize>();
#elif defined(KRIPKE_USE_HIP)
        RAJA::synchronize<RAJA::hip_synchronize>();
#endif
      }
#endif

      Kripke::Core::Set const *m_set;
      std::vector<size_t> m_chunk_to_size;
#ifndef KRIPKE_USE_CHAI
      std::vector<ElementPtr> m_chunk_to_data;
#else
      mutable std::vector<ChunkData> m_chunk_to_data;
      Kripke::ExecutionSpace m_allocation_space;
#endif
	  };

  /**
   * Defines a multi-dimensional data field defined over a Set
   */
  template<typename ELEMENT, typename ... IDX_TYPES>
  class Field : public Kripke::Core::FieldStorage<ELEMENT> {
    public:

      using Parent = Kripke::Core::FieldStorage<ELEMENT>;

      using ElementType = ELEMENT;
      static constexpr bool host_resident_normal_chai_gpu = false;
      using ElementPtr = ELEMENT*;

      static constexpr size_t NumDims = sizeof...(IDX_TYPES);

      using DefaultLayoutType = RAJA::TypedLayout<RAJA::Index_type, camp::tuple<IDX_TYPES...>>;

      //FIXME: Remove the internal namespace when new RAJA release is out
      using DefaultViewType = RAJA::internal::ViewBase<ElementType, ElementPtr, DefaultLayoutType>;
      using DeviceViewType = RAJA::internal::ViewBase<ElementType, ElementType *, DefaultLayoutType>;

      template<typename Order>
      Field(Kripke::Core::Set const &spanned_set, Order) :
        Parent(spanned_set)
      {
        setupLayouts<Order>(spanned_set);
      }

#ifdef KRIPKE_USE_CHAI
      template<typename Order>
      Field(Kripke::Core::Set const &spanned_set,
            Kripke::ExecutionSpace allocation_space,
            Order) :
        Parent(spanned_set, allocation_space)
      {
        setupLayouts<Order>(spanned_set);
      }
#endif

      template<typename Order>
      void setupLayouts(Kripke::Core::Set const &spanned_set) {
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
        auto ptr = Parent::getHostData(sdom_id);
        auto layout = m_chunk_to_layout[chunk_id];

        return DefaultViewType(ptr, layout);
      }


      RAJA_INLINE
      DeviceViewType getDeviceView(Kripke::SdomId sdom_id) const {

        size_t chunk_id = Parent::m_subdomain_to_chunk[*sdom_id];
        auto layout = m_chunk_to_layout[chunk_id];

#if defined(KRIPKE_USE_CHAI) && (defined(KRIPKE_USE_CUDA) || defined(KRIPKE_USE_HIP))
        KRIPKE_ASSERT(Parent::m_allocation_space == Kripke::GPU,
            "getDeviceView requires a GPU-backed field");
        auto ptr = Parent::getDeviceData(sdom_id);
#elif defined(KRIPKE_USE_CHAI)
        auto ptr = Parent::getHostData(sdom_id);
#else
        auto ptr = Parent::m_chunk_to_data[chunk_id];
#endif

        return DeviceViewType(ptr, layout);
      }


      template<typename Order>
      RAJA_INLINE
      auto getViewOrder(Kripke::SdomId sdom_id) const ->
        ViewType<Order, ElementType, ElementType *, IDX_TYPES...>
      {
        size_t chunk_id = Parent::m_subdomain_to_chunk[*sdom_id];

        using LInfo = LayoutInfo<Order, IDX_TYPES...>;
        using LType = typename LInfo::Layout;

        LType layout = RAJA::make_stride_one<LInfo::stride_one_dim>(m_chunk_to_layout[chunk_id]);

#if (defined(KRIPKE_USE_HIP) || defined(KRIPKE_USE_CUDA)) && defined(KRIPKE_USE_CHAI)
        if(Parent::m_allocation_space == Kripke::GPU){
          return ViewType<Order, ElementType, ElementType *, IDX_TYPES...>(Parent::getDeviceData(sdom_id), layout);
        }
        return ViewType<Order, ElementType, ElementType *, IDX_TYPES...>(Parent::getHostData(sdom_id), layout);
#elif defined(KRIPKE_USE_CHAI)
        return ViewType<Order, ElementType, ElementType *, IDX_TYPES...>(Parent::getHostData(sdom_id), layout);
#else
        return ViewType<Order, ElementType, ElementType *, IDX_TYPES...>(Parent::m_chunk_to_data[chunk_id], layout);
#endif
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
#else
        printf("  m_chunk_to_data(host): ");
        for(auto const &x : Parent::m_chunk_to_data){printf("%p ", (void*)x.host_ptr);}
        printf("\n");
        printf("  m_chunk_to_data(device): ");
        for(auto const &x : Parent::m_chunk_to_data){printf("%p ", (void*)x.device_ptr);}
        printf("\n");
#endif

        for(size_t chunk_id = 0;chunk_id < Parent::m_chunk_to_data.size();++ chunk_id){

          SdomId sdom_id(DomainVar::m_chunk_to_subdomain[chunk_id]);

          ElementType *ptr = Parent::getHostData(sdom_id);

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

  template<typename ELEMENT, bool HOST_RESIDENT_NORMAL_CHAI_GPU, typename ... IDX_TYPES>
  class FieldWithPolicy : public Kripke::Core::Field<ELEMENT, IDX_TYPES...> {
    public:
      using Parent = Kripke::Core::Field<ELEMENT, IDX_TYPES...>;
      using Parent::Parent;
      static constexpr bool host_resident_normal_chai_gpu =
          HOST_RESIDENT_NORMAL_CHAI_GPU;
  };

} } // namespace

#endif
