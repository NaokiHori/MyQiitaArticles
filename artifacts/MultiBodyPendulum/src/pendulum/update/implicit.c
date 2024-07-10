#if 0 == SCHEME

#include <math.h>
#include <float.h>
#include "pendulum.h"
#include "../internal.h"

static inline double sinc(
    const double x
){
  if (0. == x) {
    return 1.;
  } else {
    return sin(x) / x;
  }
}

static int gauss_elimination(
    const int_fast32_t nitems,
    double * restrict a,
    double * restrict x
){
  // NOTE: zero-division checks are omitted, assuming rank(A) is nitems
  // forward elimination
  for (int_fast32_t i = 0; i < nitems; i++) {
    for (int_fast32_t ii = i + 1; ii < nitems; ii++) {
      const double v = a[ii * nitems + i] / a[i * nitems + i];
      for (int_fast32_t j = i + 1; j < nitems; j++) {
        a[ii * nitems + j] -= v * a[i * nitems + j];
      }
      x[ii] -= v * x[i];
    }
  }
  // backward substitution
  for (int_fast32_t i = nitems - 1; i >= 0; i--) {
    for (int_fast32_t j = i + 1; j < nitems; j++) {
      x[i] -= a[i * nitems + j] * x[j];
    }
    x[i] /= a[i * nitems + i];
  }
  return 0;
}

static int solve_linear_system(
    const uint_t nitems,
    double * restrict a,
    double * restrict x,
    double * restrict b,
    double * restrict residual
){
  // solve A x = b
  // x stores the answer of the previous step,
  //   which are compared with the new solution to quantify the convergence
  gauss_elimination((int_fast32_t)nitems, a, b);
  // check convergence and update x
  *residual = 0.;
  for (uint_t i = 0; i < nitems; i++) {
    *residual += fabs(b[i] - x[i]);
    x[i] = b[i];
  }
  return 0;
}

static int get_new_values(
    const double v_old,
    const double a_old,
    const double dv,
    const double dt,
    double * restrict v_new,
    double * restrict a_new
){
  *v_new = v_old + dv;
  *a_new = a_old + 0.5 * (v_old + *v_new) * dt;
  return 0;
}

static int integrate(
    const double dt,
    p_object_t * restrict pendulum
){
  const uint_t nitems    = pendulum->nitems;
  double * restrict pos  = pendulum->pos;
  double * restrict vel  = pendulum->vel;
  double * restrict lhs  = pendulum->sys_lhs;
  double * restrict rhs  = pendulum->sys_rhs;
  double * restrict dvel = pendulum->dvel;
  // reset velocity increment
  for (uint_t i = 0; i < nitems; i++) {
    dvel[i] = 0.;
  }
  // iterate until converged
  const uint_t itermax = nitems << 1;
  for (uint_t iter = 0; iter < itermax; iter += 1) {
    // compute LHS array and RHS vector
    for (uint_t i = 0; i < nitems; i++) {
      const double vi_old = vel[i];
      const double ai_old = pos[i];
      const double di = dvel[i];
      double vi_new = 0.;
      double ai_new = 0.;
      get_new_values(vi_old, ai_old, di, dt, &vi_new, &ai_new);
      // RHS, potential energy contribution
      rhs[i] = (nitems - i)
        * sinc(0.5 * ai_new - 0.5 * ai_old)
        * cos(0.5 * ai_new + 0.5 * ai_old)
        * dt;
      // interactive effects
      for (uint_t j = 0; j < nitems; j++) {
        const double m = 1. * nitems - fmax(i, j);
        const double vj_old = vel[j];
        const double aj_old = pos[j];
        const double dj = dvel[j];
        double vj_new = 0.;
        double aj_new = 0.;
        get_new_values(vj_old, aj_old, dj, dt, &vj_new, &aj_new);
        // LHS array
        const double cij_old = i == j ? 1. : cos(ai_old - aj_old);
        const double cij_new = i == j ? 1. : cos(ai_new - aj_new);
        const double numer = DBL_EPSILON + (0.5 * vi_new * cij_new + 0.5 * vi_old * cij_old);
        const double denom = DBL_EPSILON + (0.5 * vi_new + 0.5 * vi_old) * (0.5 * cij_new + 0.5 * cij_old);
        const double cor = 0.5 + 0.5 * numer / denom;
        lhs[i * nitems + j] = m * cor * (0.5 * cij_new + 0.5 * cij_old);
        // RHS vector, kinetic energy contribution
        const double vj = 0.5 * vj_new + 0.5 * vj_old;
        rhs[i] -= m * vj * vj
          * sinc(0.5 * ai_new - 0.5 * ai_old - 0.5 * aj_new + 0.5 * aj_old)
          * sin(0.5 * ai_new + 0.5 * ai_old - 0.5 * aj_new - 0.5 * aj_old)
          * dt;
      }
    }
    // solve linear system A x = b and check convergence
    double residual = DBL_MAX;
    solve_linear_system(nitems, lhs, dvel, rhs, &residual);
    if (residual < DBL_EPSILON) {
      // update velocities and positions
      for (uint_t i = 0; i < nitems; i++) {
        get_new_values(vel[i], pos[i], dvel[i], dt, vel + i, pos + i);
      }
      // converged, report success
      return 0;
    }
  }
  // failed to converge with the current dt
  return 1;
}

int NSPC(update)(
    p_object_t * pendulum,
    double * restrict dt_
){
  // integrate the equations in time (from t to t+dt)
  double * restrict dt = &pendulum->dt;
  // adaptive time step size, initially try doubled value
  for (*dt = *dt * 2.; ; ) {
    if (0 == integrate(*dt, pendulum)) {
      break;
    }
    // failed to converge, retry with smaller time step size
    (*dt) *= 0.5;
    if ((*dt) < DBL_EPSILON) {
      return 1;
    }
  }
  // completed, tell how much time it proceeds to the user
  *dt_ = *dt;
  return 0;
}

#else

extern char dummy;

#endif
