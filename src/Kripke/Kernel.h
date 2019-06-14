//
// Copyright (c) 2014-19, Lawrence Livermore National Security, LLC
// and Kripke project contributors. See the COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
//

#ifndef KRIPKE_KERNEL_H__
#define KRIPKE_KERNEL_H__

#include <Kripke.h>
#include <Kripke/Core/DataStore.h>
#include <utility>

namespace Kripke {

  namespace Kernel {

    void LPlusTimes(Kripke::Core::DataStore &data_store);


    void LTimes(Kripke::Core::DataStore &data_store);


    double population(Kripke::Core::DataStore &data_store);


    void scattering(Kripke::Core::DataStore &data_store);


    void source(Kripke::Core::DataStore &data_store);


    void sweepSubdomain(Kripke::Core::DataStore &data_store, Kripke::SdomId sdom_id);


    template<typename FieldType>
    RAJA_INLINE
    void kConst(FieldType &field, Kripke::SdomId sdom_id, typename FieldType::ElementType value){
      auto view1d = field.getView1d(sdom_id);
      int num_elem = field.size(sdom_id);
      RAJA::forall<RAJA::loop_exec>(
        RAJA::RangeSegment(0, num_elem),
        [=](RAJA::Index_type i){
			 	  view1d(i) = value;
      });
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
      auto view_src = field_src.getView1d(sdom_id_src);
      auto view_dst = field_dst.getView1d(sdom_id_dst);
      int num_elem = field_src.size(sdom_id_src);

      RAJA::forall<RAJA::loop_exec>(
        RAJA::RangeSegment(0, num_elem),
        [=](RAJA::Index_type i){
          view_src(i) = view_dst(i);
      });
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
