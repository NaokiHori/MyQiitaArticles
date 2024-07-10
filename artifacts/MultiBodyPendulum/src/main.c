#include <stdio.h>
#include "pendulum.h"
#include "logger.h"
#include "visualiser.h"

int main(
    void
){
  // number of mass points
  const size_t nitems = 8;
  // max simulation time
  const double timemax = 1.e+2;
  // logging frequency
  const double log_rate = 5.e-1;
#if defined(ENABLE_XWINDOW)
  // visualiser frequency
  const double visu_rate = 1.e-1;
#endif
  // initialise state of a pendulum
  p_object_t * p_object = NULL;
  if (0 != p_solver.init(nitems, &p_object)) {
    return 1;
  }
  // initialise logger
  if (0 != logger.init(log_rate)) {
    return 1;
  }
#if defined(ENABLE_XWINDOW)
  // initialise visualiser
  if (0 != visualiser.init(visu_rate, p_object)) {
    return 1;
  }
#endif
  // main loop
  for (double time = 0.; time < timemax; ) {
    // proceed for one step
    // NOTE: dt is decided by the solver
    double dt = 0.;
    if (0 != p_solver.update(p_object, &dt)) {
      puts("failed to update pendulum");
      break;
    }
    time += dt;
    // post-process, when necessary
    if (time > timemax) {
      break;
    }
    if (time > logger.get_next_time()) {
      logger.execute(time, dt, p_object);
    }
#if defined(ENABLE_XWINDOW)
    if (time > visualiser.get_next_time()) {
      visualiser.execute(time, p_object);
    }
#endif
  }
  // clean-up objects
  p_solver.finalise(p_object);
  logger.finalise();
#if defined(ENABLE_XWINDOW)
  visualiser.finalise();
#endif
  return 0;
}

