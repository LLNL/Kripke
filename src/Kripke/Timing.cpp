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

#include <Kripke/Timing.h>

#include <Kripke.h>
#include <Kripke/Core/Comm.h>

#include <stdio.h>
#include <algorithm>
#include <vector>

using namespace Kripke;

Timing::~Timing(){
  stopAll();
}

void Timing::start(std::string const &name){
  // get or create timer
  Timer &timer = timers[name];

  // Start it up
  timer.start(name);
}

void Timing::stop(std::string const &name){
  // get or create timer
  Timer &timer = timers[name];

  // stop it
  timer.stop(name);
}

void Timing::stopAll(void){
  for(auto i : timers){
    stop(i.first);
  }
}

void Timing::print(void) const {
  Kripke::Core::Comm default_comm;
  if(default_comm.rank() != 0){
    return;
  }

  // build a sorted list of names
  std::vector<std::string> names;
  for(auto i : timers){
    names.push_back(i.first);
  }
  std::sort(names.begin(), names.end());

  std::vector<Timer const *> ord_timers;
  for(auto name : names){
    ord_timers.push_back(&timers.find(name)->second);
  }

  // Display column names
  printf("\nTimers\n");
  printf("======\n\n");
  printf("  %-16s  %12s  %12s", "Timer", "Count", "Seconds");
  printf("\n");
  printf("  ----------------  ------------  ------------\n");

  // Dislpay timer results
  for(size_t i = 0;i < names.size();++ i){
    printf("  %-16s  %12d  %12.5lf\n",
        names[i].c_str(),
        (int)ord_timers[i]->getCount(),
        ord_timers[i]->getElapsed());
  }
  
  // Now display timers in machine readable format
  printf("\n");
  printf("TIMER_NAMES:");
  for(size_t i = 0;i < names.size();++ i){
    if(i > 0){
      printf(",");
    }
    printf("%s", names[i].c_str());
  }
  printf("\n");
  printf("TIMER_DATA:");
  for(size_t i = 0;i < names.size();++ i){
    if(i > 0){
      printf(",");
    }
    printf("%lf", ord_timers[i]->getElapsed());
  }
  printf("\n");
}


double Timing::getTotal(std::string const &name) const{
  auto i = timers.find(name);
  if(i == timers.end()){
    return 0.0;
  }
  return i->second.getElapsed();
}

size_t Timing::getCount(std::string const &name) const{
  auto i = timers.find(name);
  if(i == timers.end()){
    return 0;
  }
  return i->second.getCount();
}

