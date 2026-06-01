//
// Copyright (c) 2014-25, Lawrence Livermore National Security, LLC
// and Kripke project contributors. See the Kripke/COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
//

#ifndef KRIPKE_KERNEL_H__
#define KRIPKE_KERNEL_H__

#include <Kripke.h>
#include <Kripke/ArchLayout.h>
#include <Kripke/Core/DataStore.h>
#include <utility>

#ifdef KRIPKE_USE_CHAI
#include <chai/ExecutionSpaces.hpp>
#endif

namespace Kripke {

  namespace Kernel {

    namespace detail {

      template<typename FieldType>
      RAJA_INLINE
      bool fieldUsesDevice(FieldType const &field){
#if defined(KRIPKE_USE_CHAI) && (defined(KRIPKE_USE_CUDA) || defined(KRIPKE_USE_HIP))
        return field.getAllocationSpace() == chai::GPU;
#else
        (void)field;
        return false;
#endif
      }

      template<typename FieldType>
      RAJA_INLINE
      void kConstDevice(FieldType &field, Kripke::SdomId sdom_id, typename FieldType::ElementType value){
        auto ptr = field.getDeviceData(sdom_id);
        int num_elem = field.size(sdom_id);

#if defined(KRIPKE_USE_CUDA) || defined(KRIPKE_USE_HIP)
#if defined(KRIPKE_USE_CUDA)
        RAJA::forall<RAJA::cuda_exec<256>>(
#elif defined(KRIPKE_USE_HIP)
        RAJA::forall<RAJA::hip_exec<256>>(
#else
        RAJA::forall<RAJA::seq_exec>( // should never reach this
#endif
          RAJA::RangeSegment(0, num_elem),
          KRIPKE_LAMBDA (RAJA::Index_type i){
            ptr[i] = value;
        });
#else
        (void)ptr;
        (void)num_elem;
#endif  // #if defined(KRIPKE_USE_CUDA) || defined(KRIPKE_USE_HIP)
      }

      template<typename FieldType>
      RAJA_INLINE
      void kConstHost(FieldType &field, Kripke::SdomId sdom_id, typename FieldType::ElementType value){
        auto ptr = field.getHostData(sdom_id);
        int num_elem = field.size(sdom_id);

        RAJA::forall<RAJA::seq_exec>(
          RAJA::RangeSegment(0, num_elem),
          [=](RAJA::Index_type i){
            ptr[i] = value;
        });
      }

      template<typename FieldType>
      RAJA_INLINE
      void kCopyDevice(FieldType &field_dst, Kripke::SdomId sdom_id_dst,
                       FieldType &field_src, Kripke::SdomId sdom_id_src){
        auto src = field_src.getDeviceData(sdom_id_src);
        auto dst = field_dst.getDeviceData(sdom_id_dst);
        int num_elem = field_src.size(sdom_id_src);

#if defined(KRIPKE_USE_CUDA) || defined(KRIPKE_USE_HIP)
#if defined(KRIPKE_USE_CUDA)
        RAJA::forall<RAJA::cuda_exec<256>>(
#elif defined(KRIPKE_USE_HIP)
        RAJA::forall<RAJA::hip_exec<256>>(
#else
        RAJA::forall<RAJA::seq_exec>( // should never reach this
#endif
          RAJA::RangeSegment(0, num_elem),
          KRIPKE_LAMBDA (RAJA::Index_type i){
            dst[i] = src[i];
        });
#else
        (void)src;
        (void)dst;
        (void)num_elem;
#endif  // #if defined(KRIPKE_USE_CUDA) || defined(KRIPKE_USE_HIP)
      }

      template<typename FieldType>
      RAJA_INLINE
      void kCopyHost(FieldType &field_dst, Kripke::SdomId sdom_id_dst,
                     FieldType &field_src, Kripke::SdomId sdom_id_src){
        auto src = field_src.getHostDataConst(sdom_id_src);
        auto dst = field_dst.getHostData(sdom_id_dst);
        int num_elem = field_src.size(sdom_id_src);

        RAJA::forall<RAJA::seq_exec>(
          RAJA::RangeSegment(0, num_elem),
          [=](RAJA::Index_type i){
            dst[i] = src[i];
        });
      }

    }

    void LPlusTimes(Kripke::Core::DataStore &data_store);


    void LTimes(Kripke::Core::DataStore &data_store);


    double population(Kripke::Core::DataStore &data_store);


    void scattering(Kripke::Core::DataStore &data_store);


    void source(Kripke::Core::DataStore &data_store);


    void sweepSubdomain(Kripke::Core::DataStore &data_store, Kripke::SdomId sdom_id);


    template<typename FieldType>
    RAJA_INLINE
    void kConst(FieldType &field, Kripke::SdomId sdom_id, typename FieldType::ElementType value){
      if(detail::fieldUsesDevice(field)){
        detail::kConstDevice(field, sdom_id, value);
      }
      else{
        detail::kConstHost(field, sdom_id, value);
      }
    }

    template<typename FieldType>
    RAJA_INLINE
    void kConst(FieldType &field, typename FieldType::ElementType value){
      for(Kripke::SdomId sdom_id : field.getWorkList()){
        kConst(field, sdom_id, value);
      }
    }




    template<typename FieldType>
    RAJA_INLINE
    void kCopy(FieldType &field_dst, Kripke::SdomId sdom_id_dst,
               FieldType &field_src, Kripke::SdomId sdom_id_src){
      KRIPKE_ASSERT(field_dst.size(sdom_id_dst) == field_src.size(sdom_id_src),
          "Cannot copy fields with different subdomain sizes");

      if(detail::fieldUsesDevice(field_dst) && detail::fieldUsesDevice(field_src)){
        detail::kCopyDevice(field_dst, sdom_id_dst, field_src, sdom_id_src);
      }
      else{
        detail::kCopyHost(field_dst, sdom_id_dst, field_src, sdom_id_src);
      }
    }

    template<typename FieldType>
    RAJA_INLINE
    void kCopy(FieldType &field_dst, FieldType &field_src){
      for(Kripke::SdomId sdom_id : field_dst.getWorkList()){
        kCopy(field_dst, sdom_id, field_src, sdom_id);
      }
    }

  }
}

#endif
