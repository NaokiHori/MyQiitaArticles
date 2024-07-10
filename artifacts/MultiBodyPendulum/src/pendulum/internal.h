#if !defined(PENDULUM_INTERNAL_H)
#define PENDULUM_INTERNAL_H

#include <stdint.h> // uint_fast32_t
#include "pendulum.h"

// introduce a namespace
#define NSPC(FUNC) p_solver_##FUNC

typedef uint_fast32_t uint_t;

struct p_object_t_ {
  uint_t nitems;
  // time step size, which is also a state of the system
  double dt;
  double * restrict pos;
  double * restrict vel;
  // internal buffers to avoid allocating every iteration
  double * restrict dvel;
  double * restrict sys_lhs;
  double * restrict sys_rhs;
};

extern int NSPC(init)(
    const size_t nitems_,
    p_object_t ** p_object
);

extern int NSPC(finalise)(
    p_object_t * p_object
);

extern size_t NSPC(get_nitems)(
    const p_object_t * p_object
);

extern double NSPC(get_pos)(
    const p_object_t * p_object,
    const size_t index
);

extern double NSPC(get_vel)(
    const p_object_t * p_object,
    const size_t index
);

extern int NSPC(update)(
    p_object_t * p_object,
    double * restrict dt_
);

#endif // PENDULUM_INTERNAL_H
