#if !defined(VISUALISER_H)
#define VISUALISER_H

#include "pendulum.h"

typedef struct {
  // constructor
  int (* const init)(
      const double rate,
      const p_object_t * p_object
  );
  // destructor
  int (* const finalise)(
      void
  );
  // getter
  double (* const get_next_time)(
      void
  );
  // update display
  int (* const execute)(
      const double time,
      const p_object_t * p_object
  );
} visualiser_t;

extern const visualiser_t visualiser;

#endif // VISUALISER_H
