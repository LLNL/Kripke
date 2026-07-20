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
#include <Kripke/Core/Set.h>
#include <vector>

#ifdef KRIPKE_USE_CHAI
#include <chai/ManagedArray.hpp>
#endif
#ifdef KRIPKE_USE_GPU_AWARE_MPI
#include <umpire/ResourceManager.hpp>
#include <umpire/strategy/NamedAllocationStrategy.hpp>
#endif

namespace Kripke {
namespace Core {
  template<typename ELEMENT, bool HOST_RESIDENT_NORMAL_CHAI_GPU, typename ... IDX_TYPES>
  class FieldWithPolicy;
#ifdef KRIPKE_USE_GPU_AWARE_MPI
  template<typename ELEMENT, typename ... IDX_TYPES>
  class FieldWithDirectUmpireDeviceStorage;
#endif

#ifdef KRIPKE_USE_GPU_AWARE_MPI
namespace detail {
  inline umpire::Allocator directUmpireDeviceAllocator()
  {
    auto &rm = umpire::ResourceManager::getInstance();
    char const *allocator_name = "KRIPKE_DEVICE_DIRECT";

    if(!rm.isAllocator(allocator_name)){
      return rm.makeAllocator<umpire::strategy::NamedAllocationStrategy>(
          allocator_name, rm.getAllocator("DEVICE"));
    }

    return rm.getAllocator(allocator_name);
  }
}
#endif

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
      //FIXME: Remove the internal namespace when new RAJA release is out
      using View1dType = RAJA::internal::ViewBase<ElementType, ElementPtr, Layout1dType>;


      explicit FieldStorage(Kripke::Core::Set const &spanned_set
#ifdef KRIPKE_USE_CHAI
          , chai::ExecutionSpace allocation_space = chai::CPU
#ifdef KRIPKE_USE_GPU_AWARE_MPI
          , bool direct_umpire_device_storage = false
#endif
#endif
          ) :
        m_set(&spanned_set)
#ifdef KRIPKE_USE_CHAI
        , m_allocation_space(allocation_space)
#ifdef KRIPKE_USE_GPU_AWARE_MPI
        , m_direct_umpire_device_storage(direct_umpire_device_storage && allocation_space == chai::GPU)
#endif
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

          m_chunk_to_size[chunk_id] = sdom_size;
#ifndef KRIPKE_USE_CHAI
          m_chunk_to_data[chunk_id] = new ElementType[sdom_size];
#else
#ifdef KRIPKE_USE_GPU_AWARE_MPI
          // Used only for GPU-aware MPI i/j/k_plane buffers to avoid QuickPool non-base pointers.
          if(m_direct_umpire_device_storage){
            auto &rm = umpire::ResourceManager::getInstance();
            auto host_allocator = rm.getAllocator("HOST");
            auto device_allocator = detail::directUmpireDeviceAllocator();

            m_chunk_to_data[chunk_id] = ElementPtr(
                sdom_size,
                {chai::CPU, chai::GPU}, // which CHAI execution spaces are being overridden
                {host_allocator, device_allocator}, // matching Umpire allocators for the overridden spaces {HOST, NamedAllocationStrategy}
                chai::GPU); // initial allocation space
          }
          else
#endif
          {
            m_chunk_to_data[chunk_id].allocate(sdom_size, m_allocation_space);
          }
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

#ifdef KRIPKE_USE_CHAI
        m_chunk_to_data[chunk_id].data(chai::CPU);
#endif
        ElementPtr ptr = m_chunk_to_data[chunk_id];
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
        return m_chunk_to_data[chunk_id].data(chai::CPU);
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
        return m_chunk_to_data[chunk_id].data(chai::CPU);
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
#ifdef KRIPKE_USE_GPU_AWARE_MPI
        if(m_direct_umpire_device_storage){
          return m_chunk_to_data[chunk_id].data(chai::GPU, false);
        }
#endif
#if defined(KRIPKE_USE_CUDA) || defined(KRIPKE_USE_HIP)
        return m_chunk_to_data[chunk_id].data(chai::GPU);
#else
        return m_chunk_to_data[chunk_id].data(chai::CPU);
#endif
#endif
      }

      RAJA_INLINE
      ElementType *getData(Kripke::SdomId sdom_id) const {
        return getHostData(sdom_id);
      }

#ifdef KRIPKE_USE_CHAI
      RAJA_INLINE
      chai::ExecutionSpace getAllocationSpace() const {
        return m_allocation_space;
      }

#ifdef KRIPKE_USE_GPU_AWARE_MPI
      RAJA_INLINE
      void registerDeviceTouch(Kripke::SdomId sdom_id) {
        if(m_direct_umpire_device_storage){
          return;
        }
        if(m_allocation_space == chai::GPU){
          KRIPKE_ASSERT(*sdom_id < (int)m_subdomain_to_chunk.size(),
              "sdom_id(%d) >= num_subdomains(%d)",
              (int)*sdom_id,
              (int)(int)m_subdomain_to_chunk.size());
          size_t chunk_id = m_subdomain_to_chunk[*sdom_id];
          m_chunk_to_data[chunk_id].registerTouch(chai::GPU);
        }
      }
#endif
#endif


      RAJA_INLINE
      Kripke::Core::Set const &getSet() const {
        return *m_set;
      }

    protected:
      Kripke::Core::Set const *m_set;
      std::vector<size_t> m_chunk_to_size;
      std::vector<ElementPtr> m_chunk_to_data;
#ifdef KRIPKE_USE_CHAI
      chai::ExecutionSpace m_allocation_space;
#ifdef KRIPKE_USE_GPU_AWARE_MPI
      bool m_direct_umpire_device_storage;
#endif
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

#ifdef KRIPKE_USE_GPU_AWARE_MPI
      static constexpr bool direct_umpire_device_storage = false;
#endif

#ifndef KRIPKE_USE_CHAI
      using ElementPtr = ELEMENT*;
#else
      using ElementPtr = chai::ManagedArray<ELEMENT>;
#endif

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
            chai::ExecutionSpace allocation_space,
            Order) :
        Parent(spanned_set, allocation_space)
      {
        setupLayouts<Order>(spanned_set);
      }

#ifdef KRIPKE_USE_GPU_AWARE_MPI
      template<typename Order>
      Field(Kripke::Core::Set const &spanned_set,
            chai::ExecutionSpace allocation_space,
            bool direct_umpire_device_storage_arg,
            Order) :
        Parent(spanned_set, allocation_space, direct_umpire_device_storage_arg)
      {
        setupLayouts<Order>(spanned_set);
      }
#endif
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

#ifdef KRIPKE_USE_CHAI
        Parent::m_chunk_to_data[chunk_id].data(chai::CPU);
#endif
        auto ptr = Parent::m_chunk_to_data[chunk_id];
        auto layout = m_chunk_to_layout[chunk_id];

        return DefaultViewType(ptr, layout);
      }


      RAJA_INLINE
      DeviceViewType getDeviceView(Kripke::SdomId sdom_id) const {

        size_t chunk_id = Parent::m_subdomain_to_chunk[*sdom_id];
        auto layout = m_chunk_to_layout[chunk_id];

        auto ptr = Parent::getDeviceData(sdom_id);

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
        using OrderedViewType = ViewType<Order, ElementType, ElementType *, IDX_TYPES...>;

        LType layout = RAJA::make_stride_one<LInfo::stride_one_dim>(m_chunk_to_layout[chunk_id]);

#ifndef KRIPKE_USE_CHAI
        return OrderedViewType(Parent::m_chunk_to_data[chunk_id], layout);
#else
#if defined(KRIPKE_USE_HIP) || defined(KRIPKE_USE_CUDA)
        if(Parent::m_allocation_space == chai::GPU){
          return OrderedViewType(Parent::getDeviceData(sdom_id), layout);
        }
#endif
        return OrderedViewType(Parent::m_chunk_to_data[chunk_id].data(chai::CPU), layout);
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
#ifdef KRIPKE_USE_GPU_AWARE_MPI
      static constexpr bool direct_umpire_device_storage = false;
#endif
  };

#ifdef KRIPKE_USE_GPU_AWARE_MPI
  template<typename ELEMENT, typename ... IDX_TYPES>
  class FieldWithDirectUmpireDeviceStorage : public Kripke::Core::Field<ELEMENT, IDX_TYPES...> {
    public:
      using Parent = Kripke::Core::Field<ELEMENT, IDX_TYPES...>;
      using Parent::Parent;
      static constexpr bool host_resident_normal_chai_gpu = false;
      static constexpr bool direct_umpire_device_storage = true;
  };
#endif

} } // namespace

#endif
