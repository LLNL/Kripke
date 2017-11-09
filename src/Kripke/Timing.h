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

#ifndef KRIPKE_TIMING_H__
#define KRIPKE_TIMING_H__

#include <Kripke.h>
#include <Kripke/Core/DataStore.h>

#include <RAJA/util/Timer.hpp>

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
      bool started;
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
    {
        timer.start(name);
    }
    inline ~BlockTimer(){
      timer.stop(name);
    }

    private:
      Timing &timer;
      std::string name;
  };

}

#define KRIPKE_TIMER(DS, NAME) \
  Kripke::BlockTimer BLK_TIMER_##NAME(DS.getVariable<Kripke::Timing>("timing"), #NAME);


#endif
