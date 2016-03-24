/*
 * Copyright 2010-2014 Regents of the University of Minnesota.
 * This file subject to terms of Creative Commons CC BY-NC 4.0, see LICENSE.html
 * Copyright 2014 Michael Tesch: CC BY-NC 4.0, see LICENSE.html
 *
 * author(s): Michael Tesch (tesch1@gmail.com)
 *
 */
/*! \file
 * \brief utility macros for dealing with matrices and vectors
 */

#ifndef MATVEC_H
#define MATVEC_H

#ifdef STANDARD_H
#warning "standard.h #included too early!"
#endif

#ifdef __clang__
#pragma clang diagnostic ignored "-Wunknown-pragmas"
#pragma clang diagnostic ignored "-Wgnu"
#endif

/* *******************************************************************
 * Basic data-types
 */
//#include <unistd.h>
//#include <float.h>

//#define REAL_T_IS_FLOAT
#define REAL_T_IS_DOUBLE

#ifdef REAL_T_IS_DOUBLE
#ifdef REAL_T_IS_FLOAT
#error "can't be both"
#endif
#define REAL_EPSILON DBL_EPSILON
#endif
#ifdef REAL_T_IS_FLOAT
#ifdef REAL_T_IS_DOUBLE
#error "can't be both"
#endif
#define REAL_EPSILON FLT_EPSILON
#endif

/* *
 *
 */

#define FP_EPSILON (1e-12)
#define FLT_LT(A,B) (((A) < (B)) && (fabs((A) - (B)) > FP_EPSILON))
#define FLT_GT(A,B) (((A) > (B)) && (fabs((A) - (B)) > FP_EPSILON))
#define FLT_EQ(A,B) (fabs((A) - (B)) <= FP_EPSILON)
#define FLT_NE(a, b) ((a) - (b) > FP_EPSILON || (a) - (b) < -FP_EPSILON)
#define VEC2_EQ(a, b) (FLT_EQ(a.x, b.x) && FLT_EQ(a.y, b.y))
#define VEC3_EQ(a, b) (FLT_EQ(a.x, b.x) && FLT_EQ(a.y, b.y) && FLT_EQ(a.z, b.z))

//#define FLT_NE(a, b) ((a) - (b) > FLT_EPSILON || (a) - (b) < -FLT_EPSILON)

/* *******************************************************************
 * Compound data-types
 */

/*! \brief      3x3 matrix
 */
#if defined(__OPENCL_VERSION__) && __OPENCL_VERSION__ == CL_VERSION_1_1
#ifdef REAL_T_IS_DOUBLE
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
//#pragma OPENCL EXTENSION cl_amd_fp64 : enable // subset of cl_khr_fp64
#endif
#endif


#ifndef __OPENCL_VERSION__

#if defined(HAVE_LIBOPENCL)

#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/opencl.h>
#endif

typedef cl_float2 float2;
typedef cl_float3 float3;
typedef cl_float4 float4;
typedef cl_double2 double2;
typedef cl_double3 double3;
typedef cl_double4 double4;

#else // !HAVE_LIBOPENCL

/*! \brief      1x2 vector
 */
typedef union _float2 {
  float s[2];
  struct {float x, y;};
} __attribute__((packed)) float2;

/*! \brief      1x2 vector
 */
typedef union _double2 {
  double s[2];
  struct {double x, y;};
} __attribute__((packed)) __attribute__((aligned (8))) double2;

/*! \brief      1x3 vector
 */
typedef union _float3 {
  float s[3];
  struct {float x, y, z, ww;};
} __attribute__((packed)) float3;

/*! \brief      1x3 vector
 */
typedef union _double3 {
  double s[3];
  struct {double x, y, z, ww;};
} __attribute__((packed)) __attribute__((aligned (8))) double3;

/*! \brief      1x3 vector
 */
typedef union _float4 {
  float s[4];
  union {
    float3 xyz;
    struct {float x, y, z, w;};
  };
} __attribute__((packed)) __attribute__((aligned (8))) float4;

/*! \brief      1x3 vector
 */
typedef union _double4 {
  double s[4];
  union {
    double3 xyz;
    struct {double x, y, z, w;};
  };
} __attribute__((packed)) __attribute__((aligned (8))) double4;

#endif // HAVE_LIBOPENCL
#endif // __OPENCL_VERSION__


#ifdef REAL_T_IS_FLOAT
typedef float real_t;
typedef float2 vec2_t;
typedef float3 vec3_t;
typedef float4 vec4_t;
#else
typedef double real_t;
typedef double2 vec2_t;
typedef double3 vec3_t;
typedef double4 vec4_t;
#endif

typedef struct _mat22 {
  vec2_t r1;
  vec2_t r2;
} mat22_t;

typedef struct _mat33 {
  vec3_t r1;
  vec3_t r2;
  vec3_t r3;
} mat33_t;

typedef struct _mat44 {
  vec4_t r1;
  vec4_t r2;
  vec4_t r3;
  vec4_t r4;
} mat44_t;

/* *******************************************************************
 * Matrix / Vector operators
 */

#ifdef __OPENCL_VERSION__

#define vec2_add(a,b,c) (a = (b) + (c))
#define vec2_len(x)     length(x)
#define vec2_angle(x)   acos(vec2_dot(a,b) / (vec2_len(a) * vec2_len(b)))
#define vec2_dot(x,y)   dot(x,y)
#define vec2_scale(a, s) do { (a) *= s; } while (0)

#define vec3_add(a,b,c) (a = (b) + (c))
#define vec3_len(x)     length(x)
#define vec3_angle(x)   acos(vec3_dot(a,b) / (vec3_len(a) * vec3_len(b)))
#define vec3_dot(x,y)   dot(x,y)
#define vec3_scale(a, s) do { (a) *= s; } while (0)
#define vec3_distance(a, b) distance(a, b)

#else // !__OPENCL_VERSION__

/*
 * 2D
 */

/*! \brief      dot product of two vec3_t */
#define vec2_dot(a, b) ((a).x * (b).x + (a).y * (b).y)
/*! \brief      scale a vec2_t */
#define float2_scale(a, s) do { float _SS_ = (s); (a).x *= _SS_; (a).y *= _SS_; } while (0)
/*! \brief      add two vec2_t's, store in third (a=b+c) */
#define double2_scale(a, s) do { double _SS_ = (s); (a).x *= _SS_; (a).y *= _SS_; } while (0)
/*! \brief      add two vec2_t's, store in third (a=b+c) */
#define vec2_add(a,b,c) do {                      \
    a.x = b.x + c.x;                              \
    a.y = b.y + c.y;                              \
    a.z = b.z + c.z;                              \
  } while(0)
/*! \brief      length of a vec2_t */
#define float2_len(a) sqrtf(vec2_dot(a, a))
#define double2_len(a) sqrt(vec2_dot(a, a))
/*! \brief      angle between two vec2_t's */
#define float2_angle(a,b) acosf(vec2_dot(a,b) / (float2_len(a) * float2_len(b)))
#define double2_angle(a,b) acos(vec2_dot(a,b) / (double2_len(a) * double2_len(b)))
/*! \brief      normalize a vector */
#define untyped2_normalize(a, type) ({               \
      type##2 res; type len = type##2_len(a);        \
        res.x = a.x / len;                           \
        res.y = a.y / len;                           \
        res; })
#define float2_normalize(a) untyped2_normalize(a, float)
#define double2_normalize(a) untyped2_normalize(a, double)

#ifdef REAL_T_IS_FLOAT
/*! \brief      length of a vec2_t */
#define vec2_len(a) float2_len(a)
/*! \brief      angle between two vec2_t's */
#define vec2_angle(a,b) float2_len(a)
#else
/*! \brief      length of a vec2_t */
#define vec2_len(a) double2_len(a)
/*! \brief      angle between two vec2_t's */
#define vec2_angle(a,b) double2_len(a)
#endif

/*
 * 3D
 */

/*! \brief      dot product of two vec3_t */
#define vec3_dot(a, b) ((a).x * (b).x + (a).y * (b).y + (a).z * (b).z)
/*! \brief      length of a vec3_t */
#define float3_len(a) sqrtf(vec3_dot(a, a))
#define double3_len(a) sqrt(vec3_dot(a, a))
/*! \brief      cross product of two vec3_t */
#define untyped3_cross(a, b, type) ({type res;  \
      res.x = (a).y * (b).z - (a).z * (b).y;    \
      res.y = (a).z * (b).x - (a).x * (b).z;    \
      res.z = (a).x * (b).y - (a).y * (b).x;    \
      res; })
#define float3_cross(a, b) untyped3_cross(a, b, float3)
#define double3_cross(a, b) untyped3_cross(a, b, double3)
/*! \brief      normalize a vector */
#define untyped3_normalize(a, type) ({               \
      type##3 res; type len = type##3_len(a);        \
        res.x = a.x / len;                           \
        res.y = a.y / len;                           \
        res.z = a.z / len;                           \
        res; })
#define float3_normalize(a) untyped3_normalize(a, float)
#define double3_normalize(a) untyped3_normalize(a, double)
/*! \brief      scale a vec3_t */
#define float3_scale(a, s) do { double _S_ = (s);         \
    (a).x *= _S_; (a).y *= _S_; (a).z *= _S_; } while (0)
/*! \brief      scale a vec3_t */
#define double3_scale(a, s) do { double _S_ = (s);         \
    (a).x *= _S_; (a).y *= _S_; (a).z *= _S_; } while (0)
/*! \brief      distance between two vec3_t */
#define vec3_distance(a, b)                                             \
  ({ double xx = (a).x - (b).x, yy = (a).y - (b).y, zz = (a).z - (b).z; \
    sqrt(xx*xx + yy*yy + zz*zz); })
/*! \brief      add two vec3_t's, store in third (a=b+c) */
#define untyped3_add(a, b, type) ({type res; \
      res.x = (a).x + (b).x; \
      res.y = (a).y + (b).y; \
      res.z = (a).z + (b).z; \
      res; })
#define float3_add(a, b) untyped3_add(a, b, float3)
#define double3_add(a, b) untyped3_add(a, b, double3)
/*! \brief      cross product of two vec3_t */
#define untyped3_sub(a, b, type) ({type res; \
      res.x = (a).x - (b).x; \
      res.y = (a).y - (b).y; \
      res.z = (a).z - (b).z; \
      res; })
#define float3_sub(a, b) untyped3_sub(a, b, float3)
#define double3_sub(a, b) untyped3_sub(a, b, double3)

/*! \brief      angle between two vec3_t's */
#define float3_angle(a,b) acosf(vec3_dot(a,b) / (float3_len(a) * float3_len(b)))
#define double3_angle(a,b) acos(vec3_dot(a,b) / (double3_len(a) * double3_len(b)))

/*! \brief      dot product of two vec3_t */
#define vec4_dot(a, b) ((a).x * (b).x + (a).y * (b).y + (a).z * (b).z + (a).w * (b).w)

#ifdef REAL_T_IS_FLOAT
/*! \brief      length of a vec3_t */
#define vec3_len(a) float3_len(a)
#define vec3_angle(a,b) float3_angle(a,b)
#define vec3_scale(a,s) float3_scale(a,s)
#define vec3_add(a,b) float3_add(a,b)
#define vec3_sub(a,b) float3_sub(a,b)
#define vec3_normalize(a) float3_normalize(a)
#define vec3_cross(a,b) float3_cross(a,b)
#else
/*! \brief      length of a vec3_t */
#define vec3_len(a) double3_len(a)
#define vec3_angle(a,b) double3_angle(a,b)
#define vec3_scale(a,s) double3_scale(a,s)
#define vec3_add(a,b) double3_add(a,b)
#define vec3_sub(a,b) double3_sub(a,b)
#define vec3_normalize(a) double3_normalize(a)
#define vec3_cross(a,b) double3_cross(a,b)
#endif

/*! \brief      swap two vec3_t's */
#define vec3_swap(a,b) do{ vec3_t tmp = a; a = b; b = tmp; } while(0)

/*
 * 4D
 */
/*! \brief      add two vec4_t's, store in third (a=b+c) */
#define vec4_add(a,b,c) do {                      \
    (a).x = (b).x + (c).x;                        \
    (a).y = (b).y + (c).y;                        \
    (a).z = (b).z + (c).z;                        \
    (a).w = (b).w + (c).w;                        \
  } while(0)

#endif // __OPENCL_VERSION__

#define vec2_printf(v) do {                     \
    printf("[%8.8f %8.8f]",                     \
           (v).x, (v).y);                       \
  } while (0)

#define vec3_printf(v) do {                     \
    printf("[%8.8f %8.8f %8.8f]",               \
           (v).x, (v).y, (v).z);                \
  } while (0)

#define mat33_printf(m) do {                    \
    printf("[%8.8f %8.8f %8.8f ;\n"             \
           " %8.8f %8.8f %8.8f ;\n"             \
           " %8.8f %8.8f %8.8f]\n",             \
           (m).r1.x, (m).r1.y, (m).r1.z,        \
           (m).r2.x, (m).r2.y, (m).r2.z,        \
           (m).r3.x, (m).r3.y, (m).r3.z);       \
  } while (0)

/*
 *
 * vin and vout should be unique
 */
/*! \brief      Matrix-vector inner product.  vout and vin must be different.
 *
 * \param[out]  vout    The vector receiving the resulting product.
 * \param[in]   m       The matrix to multiply.
 * \param[in]   vin     The vector to multiply.
 */
#define mat22vec2_inner(vout, m, vin) do {           \
    (vout).x = vec2_dot((m).r1, vin);                \
    (vout).y = vec2_dot((m).r2, vin);                \
  } while (0)

/*! \brief      Generate a rotation matrix for right-handed rotation around Z-axis.
 *
 * \param[out]  mat     The matrix to be populated with values.
 * \param[in]   th      The angle in radians to rotate.
 */
#define mat22_gen_rot(mat, th) do {                                    \
    real_t a, b;                                                       \
    a = cos(th);                                                       \
    b = sin(th);                                                       \
    mat.r1.x = a;      mat.r1.y = -b;                                  \
    mat.r2.x = b;      mat.r2.y = a;                                   \
  } while (0)

/*! \brief      Determinant of 2x2 matrix
 *
 * \param[in]   mat     The matrix.
 */
#define mat22_det(det, mat) do {                                \
    (det) = (mat).r1.x * (mat).r2.y - (mat).r1.y * (mat).r2.x;  \
  } while (0)

/*! \brief      Matrix-vector inner product.
 *
 * \param[out]  vout    The vector receiving the resulting product.
 * \param[in]   m       The matrix to multiply.
 * \param[in]   vin     The vector to multiply.
 */
#define mat33vec3_inner(vout, m, vin) do {             \
    vec3_t _vtmp3_;                                    \
    _vtmp3_.x = vec3_dot((m).r1, vin);                 \
    _vtmp3_.y = vec3_dot((m).r2, vin);                 \
    _vtmp3_.z = vec3_dot((m).r3, vin);                 \
    (vout) = _vtmp3_;                                  \
  } while (0)

/*! \brief      Generate a rotation matrix for rotation around X-axis.
 *
 * \param[out]  mat     The matrix to be populated with values.
 * \param[in]   a       The angle in radians to rotate.
 */
#define mat33_gen_rotX(mat, a) do {                                \
    mat.r1.x = 1; mat.r1.y = 0;       mat.r1.z = 0;                \
    mat.r2.x = 0; mat.r2.y = cos(a);  mat.r2.z = -sin(a);          \
    mat.r3.x = 0; mat.r3.y = sin(a);  mat.r3.z = cos(a);           \
  } while (0)

/*! \brief      Generate a rotation matrix for right-handed rotation around Y-axis.
 *
 * \param[out]  mat     The matrix to be populated with values.
 * \param[in]   a       The angle in radians to rotate.
 */
#define mat33_gen_rotY(mat, a) do {                                    \
    mat.r1.x = cos(a);  mat.r1.y = 0;  mat.r1.z = sin(a);              \
    mat.r2.x = 0;       mat.r2.y = 1;  mat.r2.z = 0;                   \
    mat.r3.x = -sin(a); mat.r3.y = 0;  mat.r3.z = cos(a);              \
  } while (0)

/*! \brief      Generate a rotation matrix for right-handed rotation around Z-axis.
 *
 * \param[out]  mat     The matrix to be populated with values.
 * \param[in]   th      The angle in radians to rotate.
 */
#define mat33_gen_rotZ(mat, th) do {                                   \
    double a, b;                                                       \
    a = cos(th);                                                       \
    b = sin(th);                                                       \
    mat.r1.x = a;      mat.r1.y = -b;         mat.r1.z = 0;            \
    mat.r2.x = b;      mat.r2.y = a;          mat.r2.z = 0;            \
    mat.r3.x = 0;      mat.r3.y = 0;          mat.r3.z = 1;            \
  } while (0)

/*! \brief      Determinant of 3x3 matrix
 *
 * \param[in]   mat     The matrix.
 */
#define mat33_det(det, mat) do {                                        \
    double a, b, c, d, e, f, g, h, i;                                   \
    a = (mat).r1.x;                                                     \
    b = (mat).r1.y;                                                     \
    c = (mat).r1.z;                                                     \
    d = (mat).r2.x;                                                     \
    e = (mat).r2.y;                                                     \
    f = (mat).r2.z;                                                     \
    g = (mat).r3.x;                                                     \
    h = (mat).r3.y;                                                     \
    i = (mat).r3.z;                                                     \
    (det) = a * e * i + b * f * g + c * d * h - c * e * g - b * d * i - a * f * h; \
  } while (0)

/*! \brief      Generate a matrix for right-handed "rotation around u of angle theta".
 *
 * \param[out]  mat     The matrix to be populated with values.
 * \param[in]   u       The vec3_t to rotate around - must be unit length.
 * \param[in]   th      The angle theta (radians) to rotate around u.
 *
 * http://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle
 */
#define mat33_gen_rotUT(mat, u, th) do {                                \
    double a, b, c;                                                     \
    a = cos(th);                                                        \
    b = 1.f - a;                                                        \
    c = sin(th);                                                        \
    (mat).r1.x = a +       b * u.x * u.x; (mat).r1.y = b * u.x * u.y - c * u.z; (mat).r1.z = b * u.x * u.z + c * u.y; \
    (mat).r2.x = b * u.x * u.y + c * u.z; (mat).r2.y = a       + b * u.y * u.y; (mat).r2.z = b * u.y * u.z - c * u.x; \
    (mat).r3.x = b * u.x * u.z - c * u.y; (mat).r3.y = b * u.y * u.z + c * u.x; (mat).r3.z = a +       b * u.z * u.z; \
  } while (0)

/*! \brief      Generate a matrix for right-handed "rotation around u of angle theta".
 *
 * \param[out]  mat     The matrix to be populated with values.
 * \param[in]   u       The vec3_t to rotate around - must be unit length.
 * \param[in]   th      The angle theta (radians) to rotate around u.
 *
 * http://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle
 */
#define mat44_gen_rotUT(mat, u, th) do {                                \
    double a, b, c;                                                     \
    a = cos(th);                                                        \
    b = 1.f - a;                                                        \
    c = sin(th);                                                        \
    (mat).r1.x = a +       b * u.x * u.x; (mat).r1.y = b * u.x * u.y - c * u.z; (mat).r1.z = b * u.x * u.z + c * u.y; (mat).r1.w = 0; \
    (mat).r2.x = b * u.x * u.y + c * u.z; (mat).r2.y = a       + b * u.y * u.y; (mat).r2.z = b * u.y * u.z - c * u.x; (mat).r2.w = 0; \
    (mat).r3.x = b * u.x * u.z - c * u.y; (mat).r3.y = b * u.y * u.z + c * u.x; (mat).r3.z = a +       b * u.z * u.z; (mat).r3.w = 0; \
    (mat).r4.x = (mat).r4.y = (mat).r4.z = 0; (mat).r4.w = 1;           \
  } while (0)

/*! \brief      Generate a matrix for psi,phi,theta. (ro=y,pe=x,ss=z)
 *
 * \param[out]  mat     The matrix to be populated with values.
 * \param[in]   phi     The angle phi (radians).
 * \param[in]   psi     The angle psi (radians).
 * \param[in]   theta   The angle theta (radians) (angle between z-Z).
 *
 * http://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle
 */
#define mat33_gen_euler1(mat, psi, phi, theta) do {                     \
    double a, b, c, d, e, f;                                            \
    a = sin(phi);                                                       \
    b = sin(psi);                                                       \
    c = sin(theta);                                                     \
    d = cos(phi);                                                       \
    e = cos(psi);                                                       \
    f = cos(theta);                                                     \
    (mat).r1.x = d * e; (mat).r1.y = -f * b + c * a * e; (mat).r1.z = c * b + f * a * e; \
    (mat).r2.x = d * b; (mat).r2.y = f * e + c * a * b;  (mat).r2.z = -c * e + f * a * b; \
    (mat).r3.x = -a;    (mat).r3.y = c * d;              (mat).r3.z = f * d; \
  } while (0)

/*trans
//ro : 1.000000 0.000000 0.000000
//pe : 0.000000 1.000000 0.000000
//ss : 0.000000 0.000000 1.000000
//cor
//ro : 1.000000 0.000000 0.000000
//pe : 0.000000 0.000000 1.000000
//ss : 0.000000 -1.000000 0.000000
//sag
//ro : 0.000000 1.000000 0.000000
//pe : 0.000000 0.000000 1.000000
//ss : 1.000000 0.000000 0.000000
*/

/*! \brief      Generate a matrix for psi,phi,theta. (ro=x,pe=y,ss=z)
 *
 * \param[out]  mat     The matrix to be populated with values.
 * \param[in]   phi     The angle phi (radians).
 * \param[in]   psi     The angle psi (radians).
 * \param[in]   theta   The angle theta (radians) (angle between z-Z).
 *
 * http://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle
 */
#define mat33_gen_euler2(mat, psi, phi, theta) do {                     \
    double a, b, c, d, e, f;                                            \
    a = sin(phi);                                                       \
    b = sin(psi);                                                       \
    c = sin(theta);                                                     \
    d = cos(phi);                                                       \
    e = cos(psi);                                                       \
    f = cos(theta);                                                     \
    (mat).r1.x = a * e - d * f * b;  (mat).r1.y = -a * b - d * f * e; (mat).r1.z = c * d; \
    (mat).r2.x = -d * e - a * f * e; (mat).r2.y = d * b - a * f * e;  (mat).r2.z = c * a; \
    (mat).r3.x = b * c;              (mat).r3.y = e * c;              (mat).r3.z = f; \
    mat33_t tmp = mat;                                                  \
    mat33_transpose(mat, tmp);                                          \
  } while (0)

/*trans
//ro : 0.000000 -1.000000 0.000000
//pe : -1.000000 0.000000 0.000000
//ss : 0.000000 0.000000 1.000000
//cor
//ro : 0.000000 0.000000 1.000000
//pe : -1.000000 0.000000 0.000000
//ss : 0.000000 1.000000 0.000000
//sag
//ro : 0.000000 0.000000 1.000000
//pe : 0.000000 1.000000 0.000000
//ss : 1.000000 0.000000 0.000000
*/

#define mat33_transpose(mat, ori) do {                                  \
    (mat).r1.x = (ori).r1.x;  (mat).r1.y = (ori).r2.x;  (mat).r1.z = (ori).r3.x; \
    (mat).r2.x = (ori).r1.y;  (mat).r2.y = (ori).r2.y;  (mat).r2.z = (ori).r3.y; \
    (mat).r3.x = (ori).r1.z;  (mat).r3.y = (ori).r2.z;  (mat).r3.z = (ori).r3.z; \
  } while (0)

#define mat33_invert(matINV, mat) do {                                  \
    double det;                                                         \
    mat33_det(det, mat);                                                \
    if (FLT_NE(det, 0.))                                                \
      det = 1. / det;                                                   \
    else                                                                \
      det = 0.;                                                         \
    double A, B, C, D, E, F, G, H, I;                                   \
    double                                                              \
      a, b, c,                                                          \
      d, e, f,                                                          \
      g, h, i;                                                          \
    a = (mat).r1.x;                                                     \
    b = (mat).r1.y;                                                     \
    c = (mat).r1.z;                                                     \
    d = (mat).r2.x;                                                     \
    e = (mat).r2.y;                                                     \
    f = (mat).r2.z;                                                     \
    g = (mat).r3.x;                                                     \
    h = (mat).r3.y;                                                     \
    i = (mat).r3.z;                                                     \
    A = e * i - f * h;                                                  \
    B = f * g - d * i;                                                  \
    C = d * h - e * g;                                                  \
    D = c * h - b * i;                                                  \
    E = a * i - c * g;                                                  \
    F = g * b - a * h;                                                  \
    G = b * f - c * e;                                                  \
    H = c * d - a * f;                                                  \
    I = a * e - b * d;                                                  \
    (matINV).r1.x = det * A;  (matINV).r1.y = det * D;  (matINV).r1.z = det * G; \
    (matINV).r2.x = det * B;  (matINV).r2.y = det * E;  (matINV).r2.z = det * H; \
    (matINV).r3.x = det * C;  (matINV).r3.y = det * F;  (matINV).r3.z = det * I; \
  } while (0)

#define mat33_eq(eq, mat1, mat2) do {            \
    if (FLT_NE((mat1).r1.x, (mat2).r1.x) ||      \
        FLT_NE((mat1).r1.y, (mat2).r1.y) ||      \
        FLT_NE((mat1).r1.z, (mat2).r1.z) ||      \
        FLT_NE((mat1).r2.x, (mat2).r2.x) ||      \
        FLT_NE((mat1).r2.y, (mat2).r2.y) ||      \
        FLT_NE((mat1).r2.z, (mat2).r2.z) ||      \
        FLT_NE((mat1).r3.x, (mat2).r3.x) ||      \
        FLT_NE((mat1).r3.y, (mat2).r3.y) ||      \
        FLT_NE((mat1).r3.z, (mat2).r3.z))        \
      eq = 0;                                    \
    else                                         \
      eq = 1;                                    \
  } while (0)
        
#define mat33_clip1(v) if (fabs(v) < 1e-14) v = 0.0;
#define mat33_clip(mat) do {                                            \
    mat33_clip1((mat).r1.x);  mat33_clip1((mat).r1.y);  mat33_clip1((mat).r1.z); \
    mat33_clip1((mat).r2.x);  mat33_clip1((mat).r2.y);  mat33_clip1((mat).r2.z); \
    mat33_clip1((mat).r3.x);  mat33_clip1((mat).r3.y);  mat33_clip1((mat).r3.z); \
  } while (0)

#define mat33_zero(mat) do {                                          \
    bzero((mat).m, sizeof(mat));                                      \
  } while (0)

#endif // MATVEC_H

