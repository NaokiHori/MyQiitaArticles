#include "pendulum.h"
#include "internal.h"

const p_solver_t p_solver = {
  .init       = NSPC(init),
  .finalise   = NSPC(finalise),
  .get_nitems = NSPC(get_nitems),
  .get_pos    = NSPC(get_pos),
  .get_vel    = NSPC(get_vel),
  .update     = NSPC(update),
};

