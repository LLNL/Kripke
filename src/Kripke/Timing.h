//
// Copyright (c) 2014-19, Lawrence Livermore National Security, LLC
// and Kripke project contributors. See the COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
//

#ifndef KRIPKE_TIMING_H__
#define KRIPKE_TIMING_H__

#include <Kripke.h>
#include <Kripke/Core/DataStore.h>

#include <RAJA/util/Timer.hpp>

#ifdef KRIPKE_USE_CALIPER
#include <caliper/cali.h>
#endif

#include <string>
#include <map>

namespace Kripke {

  class Timer {
    public:
      RAJA_INLINE
      Timer() :
        started(false),
        elapsed(0.),
        count(0)
      {}

      RAJA_INLINE
      void start(std::string const &my_name) {
        timer.stop(my_name.c_str());
        timer.reset();
        timer.start(my_name.c_str());
        started = true;
        ++ count;
      }

      RAJA_INLINE
      void stop(std::string const &my_name) {
        if(started){
          timer.stop(my_name.c_str());
          elapsed += timer.elapsed();
        }
      }

      RAJA_INLINE
      size_t getCount() const {
        return count;
      }

      RAJA_INLINE
      double getElapsed() const {
        return elapsed;
      }

    private:
      bool   started;
      double elapsed;
      size_t count;
      RAJA::Timer timer;
  };

  class Timing : public Kripke::Core::BaseVar {
    public:
      virtual ~Timing();

      void start(std::string const &name);
      void stop(std::string const &name);

      void stopAll(void);

      void print(void) const;
      double getTotal(std::string const &name) const;
      size_t getCount(std::string const &name) const;

    private:
      using TimerMap = std::map<std::string, Timer>;
      TimerMap timers;
  };


  // Aides timing a block of code, with automatic timer stopping
  class BlockTimer {
    public:
    inline BlockTimer(Timing &timer_obj, std::string const &timer_name) :
        timer(timer_obj),
        name(timer_name)
#ifdef KRIPKE_USE_CALIPER
        , cali_timer(cali_annot.begin(timer_name.c_str()))
#endif
    {
        timer.start(name);
    }
    inline ~BlockTimer(){
      timer.stop(name);
    }

    private:
      Timing &timer;
      std::string name;
#ifdef KRIPKE_USE_CALIPER
      cali::Annotation::Guard cali_timer;
      static cali::Annotation cali_annot;
#endif
  };

}

#define KRIPKE_TIMER(DS, NAME) \
  Kripke::BlockTimer BLK_TIMER_##NAME(DS.getVariable<Kripke::Timing>("timing"), #NAME);

#endif
