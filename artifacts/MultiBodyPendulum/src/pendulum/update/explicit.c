#if 1 == SCHEME

#include <math.h>
#include <float.h>
#include "pendulum.h"
#include "../internal.h"

static int cholesky_decomposition(
    const int_fast32_t nitems,
    double * restrict a,
    double * restrict x
){
  // solve A x = b by means of the Cholesky decomposition (LDLT decomposition),
  //   where x stores the input and output
  // L D LT x = b
  //   D LT x = L^{-1}  b = y
  //     LT x = D^{-1}  y = z
  //        x = LT^{-1} z
  // NOTE: zero-division checks are omitted, assuming rank(A) is nitems
  // LDLT decomposition
  for (int_fast32_t i = 0; i < nitems; i++) {
    // lower triangular matrix
    for (int_fast32_t j = 0; j < i; j++) {
      for (int_fast32_t k = 0; k < j; k++) {
        const double lik = a[i * nitems + k];
        const double ljk = a[j * nitems + k];
        const double dk  = a[k * nitems + k];
        a[i * nitems + j] -= lik * ljk * dk;
      }
      a[i * nitems + j] /= a[j * nitems + j];
    }
    // diagonal component
    for (int_fast32_t j = 0; j < i; j++) {
      const double dj  = a[j * nitems + j];
      const double lij = a[i * nitems + j];
      a[i * nitems + i] -= lij * lij * dj;
    }
  }
  // eliminations
  // L y = b
  for (int_fast32_t i = 0; i < nitems; i++) {
    for (int_fast32_t j = 0; j < i; j++) {
      x[i] -= a[i * nitems + j] * x[j];
    }
  }
  // D z = y
  for (int_fast32_t i = 0; i < nitems; i++) {
    x[i] /= a[i * nitems + i];
  }
  // LT x = z
  for (int_fast32_t i = nitems - 1; i >= 0; i--) {
    for (int_fast32_t j = i + 1; j < nitems; j++) {
      x[i] -= a[j * nitems + i] * x[j];
    }
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
  cholesky_decomposition((int_fast32_t)nitems, a, b);
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
  *a_new = a_old + *v_new * dt;
  return 0;
}

static int integrate(
    const double dt,
    p_object_t * pendulum
){
  const uint_t nitems    = pendulum->nitems;
  double * restrict pos  = pendulum->pos;
  double * restrict vel  = pendulum->vel;
  double * restrict lhs  = pendulum->sys_lhs;
  double * restrict rhs  = pendulum->sys_rhs;
  double * restrict dvel = pendulum->dvel;
  // compute the LHS array
  for (uint_t i = 0; i < nitems; i++) {
    for (uint_t j = 0; j < nitems; j++) {
      const double m = 1. * nitems - fmax(i, j);
      lhs[i * nitems + j] = m * cos(pos[i] - pos[j]);
    }
  }
  // compute the RHS vector
  for (uint_t i = 0; i < nitems; i++) {
    // potential energy contribution
    rhs[i] = (nitems - i) * cos(pos[i]) * dt;
    // kinetic energy contribution
    for (uint_t j = 0; j < nitems; j++) {
      const double m = 1. * nitems - fmax(i, j);
      rhs[i] -= m * vel[j] * vel[j] * sin(pos[i] - pos[j]) * dt;
    }
  }
  // solve linear system A x = b
  double residual = DBL_MAX;
  solve_linear_system(nitems, lhs, dvel, rhs, &residual);
  // update velocities and positions
  for (uint_t i = 0; i < nitems; i++) {
    get_new_values(vel[i], pos[i], dvel[i], dt, vel + i, pos + i);
  }
  return 0;
}

int NSPC(update)(
    p_object_t * pendulum,
    double * restrict dt_
){
  // integrate the equations in time (from t to t+dt)
  double * restrict dt = &pendulum->dt;
  // NOTE: use ad-hoc fixed time step size
  *dt = 1.e-3 / pendulum->nitems;
  if (0 != integrate(*dt, pendulum)) {
    return 1;
  }
  // completed, tell how much time it proceeds to the user
  *dt_ = *dt;
  return 0;
}

#else

extern char dummy;

#endif
