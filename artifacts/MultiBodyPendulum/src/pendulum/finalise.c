#include "memory.h"
#include "pendulum.h"
#include "internal.h"

int NSPC(finalise)(
    p_object_t * p_object
){
  memory_free(p_object->pos);
  memory_free(p_object->vel);
  memory_free(p_object->sys_lhs);
  memory_free(p_object->sys_rhs);
  memory_free(p_object->dvel);
  memory_free(p_object);
  return 0;
}

