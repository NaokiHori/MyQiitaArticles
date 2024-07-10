#if !defined(PENDULUM_H)
#define PENDULUM_H

#include <stddef.h> // size_t

// serves as the opaque pointer to the state of a pendulum
typedef struct p_object_t_ p_object_t;

// methods to manipulate p_object
typedef struct {
  // constructor
  int (* const init) (
      const size_t nitems,
      p_object_t ** pendulum
  );
  // destructor
  int (* const finalise) (
      p_object_t * pendulum
  );
  // getters
  // get number of objects
  size_t (* const get_nitems) (
      const p_object_t * pendulum
  );
  // get position (theta) of n-th object
  double (* const get_pos) (
      const p_object_t * pendulum,
      const size_t index
  );
  // get velocity (omega) of n-th object
  double (* const get_vel) (
      const p_object_t * pendulum,
      const size_t index
  );
  // update velocities and positions
  int (* const update) (
      p_object_t * pendulum,
      double * dt
  );
} p_solver_t;

extern const p_solver_t p_solver;

#endif // PENDULUM_H
