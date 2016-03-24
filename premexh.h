/* this file exists because mex.h doesn't always play well
 * with the host system, and needs some stuff defined before 
 * it can be included.
 */
#if !defined(__APPLE__) && 1
#include <uchar.h>
#else
typedef signed char char16_t;
#endif
