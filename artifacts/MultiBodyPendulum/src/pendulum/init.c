#include <stdio.h>
#include <math.h>
#include "memory.h"
#include "pendulum.h"
#include "internal.h"

int NSPC(init)(
    const size_t nitems_,
    p_object_t ** p_object
){
  const size_t nitems_max = 1 << 8;
  if(nitems_ >= nitems_max){
    printf("too big nitems: %zu, which should be smaller than %zu\n", nitems_, nitems_max);
    return 1;
  }
  const uint_t nitems = (uint_t)(nitems_);
  (*p_object) = memory_calloc(1, sizeof(p_object_t));
  // assign scalars
  // NOTE: an arbitrary positive number can be used for dt
  (*p_object)->nitems = (uint_t)nitems_;
  (*p_object)->dt = 1.;
  // allocate (and assign if needed) vectors
  double * restrict * pos     = &(*p_object)->pos;
  double * restrict * vel     = &(*p_object)->vel;
  double * restrict * sys_lhs = &(*p_object)->sys_lhs;
  double * restrict * sys_rhs = &(*p_object)->sys_rhs;
  double * restrict * dvel    = &(*p_object)->dvel;
  *pos     = memory_calloc(         nitems, sizeof(double));
  *vel     = memory_calloc(         nitems, sizeof(double));
  *sys_lhs = memory_calloc(nitems * nitems, sizeof(double));
  *sys_rhs = memory_calloc(         nitems, sizeof(double));
  *dvel    = memory_calloc(         nitems, sizeof(double));
  // give initial condition
  // give velocity such that
  //   vertically-aligned state is one (unstable) solution
  const double v = sqrt(6. / (2 * nitems + 1));
  for(size_t n = 0; n < nitems; n++){
    (*pos)[n] = 0.;
    (*vel)[n] = v;
  }
  return 0;
}

