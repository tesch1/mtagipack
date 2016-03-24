/*
 * Copyright 2010-2014 Regents of the University of Minnesota.
 * This file subject to terms of Creative Commons CC BY-NC 4.0, see LICENSE.html
 *
 * author(s): Michael Tesch (tesch1@gmail.com)
 *
 */
/*! \file       matvec tests
 * \brief       test utility macros for dealing with matrices and vectors
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matvec.h"

#define xTEST(x,X) do {                         \
    if ((X) != (x))                             \
      printf("failed: " #X " %d\n", __LINE__);  \
    else                                        \
      printf(" ok: " #X "\n");                  \
  } while (0)

int main(int argc, char *argv[])
{
  xTEST(1,FLT_LT(-1,-.99));
  xTEST(1,FLT_GT(-1,-1.99));
  xTEST(1,FLT_EQ(-1,-1+FP_EPSILON/3));
  xTEST(1,FLT_EQ(1,1+FP_EPSILON/3));
  xTEST(1,FLT_EQ(1,1-FP_EPSILON/3));
  xTEST(1,FLT_NE(1,1+FP_EPSILON));
  xTEST(0,FLT_EQ(1,1+FP_EPSILON));

  double det;

  vec2_t x2 = {{1,0}};
  vec2_t y2 = {{0,1}};
  vec2_t res2;
  mat22_t m2;
  mat22_gen_rot(m2, M_PI/2);
  mat22_det(det, m2);
  xTEST(1,FLT_EQ(det,1));
  mat22vec2_inner(res2, m2, x2);
  xTEST(1,VEC2_EQ(y2, res2)); vec2_printf(res2); printf("\n");

  vec3_t x3 = {{1,0,0}};
  vec3_t y3 = {{0,1,0}};
  vec3_t z3 = {{0,0,1}};
  vec3_t res3;
  mat33_t m3, m3i, m3t;

  // x around z
  mat33_gen_rotUT(m3, z3, M_PI/2);
  mat33_det(det, m3);
  xTEST(1,FLT_EQ(det,1));
  mat33_invert(m3i, m3);
  mat33_transpose(m3t, m3);
  mat33_eq(det, m3i, m3t);
  xTEST(1,det); mat33_printf(m3); mat33_printf(m3i); mat33_printf(m3t);
  mat33vec3_inner(res3, m3, x3);
  xTEST(1,VEC3_EQ(y3, res3)); vec3_printf(res3); printf("\n");

  // y around x
  mat33_gen_rotUT(m3, x3, M_PI/2);
  mat33_det(det, m3);
  xTEST(1,FLT_EQ(det,1));
  mat33_invert(m3i, m3);
  mat33_transpose(m3t, m3);
  mat33_eq(det, m3i, m3t);
  xTEST(1,det);
  mat33vec3_inner(res3, m3, y3);
  xTEST(1,VEC3_EQ(z3, res3)); vec3_printf(res3); printf("\n");

  // z around -x
  mat33_gen_rotUT(m3, x3, M_PI/2);
  mat33_det(det, m3);
  xTEST(1,FLT_EQ(det,1));
  mat33_invert(m3i, m3);
  mat33_transpose(m3t, m3);
  mat33_eq(det, m3i, m3t);
  xTEST(1,det);
  mat33vec3_inner(res3, m3, z3);
  mat33vec3_inner(res3, m3, res3);
  mat33vec3_inner(res3, m3, res3);
  xTEST(1,VEC3_EQ(y3, res3)); vec3_printf(res3); printf("\n");

  // z around y
  mat33_gen_rotUT(m3, y3, M_PI/2);
  mat33_det(det, m3);
  xTEST(1,FLT_EQ(det,1));
  mat33_invert(m3i, m3);
  mat33_transpose(m3t, m3);
  mat33_eq(det, m3i, m3t);
  xTEST(1,det);
  mat33vec3_inner(res3, m3, z3);
  xTEST(1,VEC3_EQ(x3, res3)); vec3_printf(res3); printf("\n");

  int fail = 0;
  for (int ii = 0; ii < 100; ii++) {
    res3.x = drand48();
    res3.y = drand48();
    res3.z = drand48();
    vec3_scale(res3, 1.0 / vec3_len(res3));
    double theta = drand48() * 4 * M_PI;
    mat33_gen_rotUT(m3, res3, theta);
    mat33_det(det, m3);
    if (!FLT_EQ(det,1)) { fail = 1; printf("fail: det not unitary: "); vec3_printf(res3); printf(" %f\n", theta); }
    mat33_invert(m3i, m3);
    mat33_transpose(m3t, m3);
    mat33_eq(det, m3i, m3t);
    if (!FLT_EQ(det,1)) { fail = 1; printf("fail: inv != trans: "); vec3_printf(res3); printf(" %f\n", theta); }
  }
  if (!fail) printf(" ok!\n");
}
