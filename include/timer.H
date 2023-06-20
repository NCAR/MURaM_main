#ifndef __TIMER_INCLUDED__
#define __TIMER_INLUCDED__

#include <time.h>

class Timer{
private:
  time_t  t0; // Wall clock initialization time
  clock_t c0; // CPU  clock initialization time
public:
  enum ttype{CPU,WALL};
  Timer() { time(&t0); c0 = clock(); }
  void set() { time(&t0); c0 = clock(); }
  double get(const ttype tp=CPU) {
    switch (tp) {
    case WALL:
      time_t t1; time(&t1);
      return difftime(t1,t0);
    case CPU:
    default:
      return ((double)(clock()-c0))/CLOCKS_PER_SEC;
    }
  }
  double cpu() {
    return ((double)(clock()-c0))/CLOCKS_PER_SEC;
  }
  double wall() {
    time_t t1; time(&t1);
    return difftime(t1,t0);
  }
};

#endif
