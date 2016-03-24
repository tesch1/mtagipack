/*
 * Copyright 2010-2014 Regents of the University of Minnesota.
 * This file subject to terms of Creative Commons CC BY-NC 4.0, see LICENSE.html
 * Copyright 2012 Michael Tesch.  All rights reserved.
 *
 * author(s): Michael Tesch (tesch1@gmail.com),
 */
/*! \file
 * \brief	Load agilent parameter file into a struct
 *
 * \section     Building
 *
 * Use mtmatpack_setup.m or mtagipack_setup.m
 *
 * \section     Example
\verbatim
octave:2> mexLoadAgilentParams
error: mexLoadAgilentParams: usage: paramstruct = mexLoadAgilentParams(parfile)

octave:2> params = mexLoadAgilentParams('/Users/tesch/Documents/umn/data/4T/rcswift_normal_01.fid/procpar');
octave:4> params.nv
ans =  4096
octave:5> params.sw
ans =  8.3333e+05
octave:7> params.saveglobal_
ans =

{
  [1,1] = probe
  [2,1] = lcpeak
  [3,1] = loc
  [4,1] = lockpower
  [5,1] = lockgain
  [6,1] = lockphase
  [7,1] = lockfreq
  [8,1] = z0
  [9,1] = lkof
  [10,1] = vloc
  [11,1] = vrack
  [12,1] = vzone
  [13,1] = vproto
  [14,1] = pkpick
  [15,1] = parstyle
  [16,1] = operator
  [17,1] = studyid
  [18,1] = systemname
  [19,1] = probetype
}

octave:8> params.saveglobal_{8}
ans = z0
octave:9> 
octave:9> params.dx
ans =

   0
   1
   1
   0
  -1
  -1
   0

octave:10> params.dx(5)
ans = -1

\endverbatim

 */
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include "premexh.h"
#include "mex.h"
#include "cprocpar.h"

#ifdef DOXYGEN
#define mexFunction mexLoadAgilentParams
#endif

/* paramstruct = mexLoadAgilentParams(filename) */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  char measfile[512];
  char * errmsg = 
    "usage: paramstruct = mexLoadAgilentParams(parfile)\n";
  procpar_t * pp;
  int count;
  int i;
  const char ** keys;
  mwSize dims[2] = {1, 1};

  if (nrhs != 1 || nlhs != 1 || !mxIsChar(prhs[0])) {
    mexWarnMsgTxt(" bad params");
    mexErrMsgTxt(errmsg);
  }

  if (mxGetString(prhs[0], measfile, sizeof(measfile))) {
    mexWarnMsgTxt(" can't get param file name from argument");
    mexErrMsgTxt(errmsg);
  }

  /* read the procpar */
  pp = pp_open(measfile, NULL, TREE_PROC);
  if (!pp) {
    mexWarnMsgTxt(measfile);
    mexErrMsgTxt(":unable to open param file");
  }
  count = pp_size(pp);

  /* gather param names into key array */
  keys = malloc(sizeof(char *) * count);
  if (!keys)
    mexErrMsgTxt("malloc");
  for (i = 0; i < count; i++) {
    keys[i] = pp_get_next(pp);
  }

  /* Construct outdata */
  plhs[0] = mxCreateStructArray(2, dims, count, keys);

  /* set struct values */
  for (i = 0; i < count; i++) {
    mxArray *field_value = NULL; // compiler warning
    int j, fieldcount;

    fieldcount = pp_get_count(pp, keys[i]);
    dims[0] = fieldcount;

    switch (pp_get_type(pp, keys[i])) {
    case VAR_INT:
      field_value = mxCreateNumericArray(2, dims, mxINT32_CLASS, mxREAL);
      for (j = 0; j < fieldcount; j++) {
        double val;
        if (pp_get_real(pp, keys[i], j, &val))
          mexErrMsgTxt("mexLoadAgilentParams error in pp_get_real() (INT)");
        ((int *)mxGetPr(field_value))[j] = val;
      }
      break;
    case VAR_REAL:
      field_value = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
      for (j = 0; j < fieldcount; j++) {
        double val;
        if (pp_get_real(pp, keys[i], j, &val))
          mexErrMsgTxt("mexLoadAgilentParams error in pp_get_real()");
        mxGetPr(field_value)[j] = val;
        //printf("%s=%f\n", keys[i], pp_get_real(pp, keys[i], j));
      }
      break;
    case VAR_STRING:
      field_value = mxCreateCellMatrix(fieldcount, 1);
      for (j = 0; j < fieldcount; j++) {
        const char * strptr;
        if (pp_get_string(pp, keys[i], j, &strptr))
          mexErrMsgTxt("mexLoadAgilentParams error in pp_get_string()");
        mxSetCell(field_value, j, mxCreateString(strptr));
      }
      break;
    case VAR_NOEXIST:
    default:
      mexErrMsgTxt("mexLoadAgilentParams internal error");
    }
    mxSetFieldByNumber(plhs[0], 0, i, field_value);
  }
  //free(keys);
}
