#include <stdio.h>
#include <stdlib.h>
#include "pendulum.h"
#include "internal.h"

static void sanitise_index(
    const p_object_t * p_object,
    const size_t index
){
  const size_t nitems = p_object->nitems;
  if(index >= nitems){
    printf("%s: out-of-bounds, index: %zu, while nitems: %zu\n", __func__, index, nitems);
    exit(EXIT_FAILURE);
  }
}

size_t NSPC(get_nitems)(
    const p_object_t * p_object
){
  return p_object->nitems;
}

double NSPC(get_pos)(
    const p_object_t * p_object,
    const size_t index
){
  sanitise_index(p_object, index);
  return p_object->pos[index];
}

double NSPC(get_vel)(
    const p_object_t * p_object,
    const size_t index
){
  sanitise_index(p_object, index);
  return p_object->vel[index];
}

