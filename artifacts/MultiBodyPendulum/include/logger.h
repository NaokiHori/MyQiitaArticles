#if !defined(LOGGER_H)
#define LOGGER_H

#include "pendulum.h"

typedef struct {
  // constructor
  int (* const init)(
      const double rate
  );
  // destructor
  int (* const finalise)(
      void
  );
  // getter
  double (* const get_next_time)(
      void
  );
  // output values to be monitored
  int (* const execute)(
      const double time,
      const double dt,
      const p_object_t *pendulum
  );
} logger_t;

extern const logger_t logger;

#endif // LOGGER_H
