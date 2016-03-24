#ifndef DEBUG_H
#define DEBUG_H 1

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <assert.h>

#ifdef __cplusplus
#error "use debug.hpp"
#endif

#ifndef IF_DIE // sometimes deffined in debug.hxx

//! print error and exit(-1) if condition is true
#define IF_DIE(cond, msg...) do {                                       \
    if ((cond)) {                                                       \
      fprintf(stderr, "%s:%d : assert failed: '%s'\n", __FILE__, __LINE__, #cond); \
      if (errno) perror("");                                            \
      fprintf(stderr, msg);                                             \
      fprintf(stderr, ", dying\n");                                     \
      fprintf(stderr, "\n");                                            \
      assert(0);                                                        \
    }                                                                   \
  } while (0)
#endif

//! print error and return -1 if condition is true
#define IF_FAIL(cond, ret, msg...) do {                                 \
    if ((cond)) {                                                       \
      fprintf(stderr, "%s:%d : failed: '%s'\n", __FILE__, __LINE__, #cond); \
      if (errno) perror("");                                            \
      fprintf(stderr, msg);                                             \
      fprintf(stderr, ", returning\n");                                 \
      fprintf(stderr, "\n");                                            \
      return (ret);                                                     \
    }                                                                   \
  } while (0)

#endif // DEBUG_H
