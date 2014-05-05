/******************************************************************************
 *
 * Header file for doing timing
 *
 *****************************************************************************/

#ifndef KRIPKE_TIMING_H__
#define KRIPKE_TIMING_H__

#include<string>
#include<map>

struct Timer {
  Timer() :
    started(false),
    start_time(0.0),
    total_time(0.0),
    count(0)
  {}
  
  bool started;
  double start_time;
  double total_time;
  size_t count;
};

class Timing {
  public:
    void start(std::string const &name);
    void stop(std::string const &name);
    
    void stopAll(void);
    void clear(void);
    
    void print(void) const;
    
  private:
    typedef std::map<std::string, Timer> TimerMap;
    TimerMap timers;
};


#endif
