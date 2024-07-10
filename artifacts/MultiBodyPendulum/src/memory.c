#include <stdio.h>
#include <stdlib.h>
#include "memory.h"

void * memory_calloc(
    const size_t count,
    const size_t size
){
  void * ptr = calloc(count, size);
  if (NULL == ptr) {
    // fatal, abort program
    fprintf(stderr, "memory allocation error: %zu, %zu\n", count, size);
    exit(EXIT_FAILURE);
  }
  return ptr;
}

void memory_free(
    void * ptr
){
  free(ptr);
}

