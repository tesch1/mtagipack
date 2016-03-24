/*
 * Copyright 2010-2014 Regents of the University of Minnesota.
 * This file subject to terms of Creative Commons CC BY-NC 4.0, see LICENSE.html
 * Copyright 2012 Michael Tesch.  All rights reserved.
 *
 * author(s): Michael Tesch (tesch1@gmail.com),
 */
/*! \file
 * \brief	Load agilent traces into an array.
 *
 * \section     Building
 *
 * Use mtmatpack_setup.m or mtagipack_setup.m
 *
 * \section     Example
\verbatim
>> mexLoadAgilentTraces
Error using mexLoadAgilentTraces
usage: [traces mainh blockh] = mexLoadAgilentTraces(filename)

>> [t, mainh, blockh] = mexLoadAgilentTraces('/home/tesch/data/data.fid');
>> size(t)

ans =

         256        3000           8

>> mainh

mainh = 

      nblocks: 8
      ntraces: 3000
           np: 256
       ebytes: 4
       tbytes: 2048
       bbytes: 6144028
       status: 69
      vers_id: 0
    nbheaders: 1

>> blockh

blockh = 

     status: [8x1 uint16]
      scale: [8x1 uint16]
       mode: [8x1 uint16]
      index: [8x1 uint16]
    ctcount: [8x1 int32]
      lpval: [8x1 single]
      rpval: [8x1 single]
        lvl: [8x1 single]
        tlt: [8x1 single]

>> 
\endverbatim

 */
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "premexh.h"
#include "mex.h"
#include "varian.h"

#include "mexCommons.h" // for mexStructToArray()

#define NARRAY(a)       ((sizeof a) / sizeof (a)[0])
#define MH_FIELD_DEF(fname, ftype) {#fname, (size_t)&((vmain_header_t*)NULL)->fname, \
      sizeof(((vmain_header_t*)NULL)->fname), ftype}
#define BH_FIELD_DEF(fname, ftype) {#fname, (size_t)&((vblock_header_t*)NULL)->fname, \
      sizeof(((vblock_header_t*)NULL)->fname), ftype}

/*! descriptions of all fields in an varian main header
 */
static struct_fielddesc_t mheader_fields[] = {
  MH_FIELD_DEF(nblocks, TY_S32),
  MH_FIELD_DEF(ntraces, TY_S32),
  MH_FIELD_DEF(np, TY_S32),
  MH_FIELD_DEF(ebytes, TY_S32),
  MH_FIELD_DEF(tbytes, TY_S32),
  MH_FIELD_DEF(bbytes, TY_S32),
  MH_FIELD_DEF(status, TY_U16),
  MH_FIELD_DEF(vers_id, TY_U16),
  MH_FIELD_DEF(nbheaders, TY_S32),
};

/*! descriptions of all fields in an varian block header
 */
static struct_fielddesc_t bheader_fields[] = {
  BH_FIELD_DEF(status, TY_U16),
  BH_FIELD_DEF(scale, TY_U16),
  BH_FIELD_DEF(mode, TY_U16),
  BH_FIELD_DEF(index, TY_U16),
  BH_FIELD_DEF(ctcount, TY_S32),
  BH_FIELD_DEF(lpval, TY_FLOAT),
  BH_FIELD_DEF(rpval, TY_FLOAT),
  BH_FIELD_DEF(lvl, TY_FLOAT),
  BH_FIELD_DEF(tlt, TY_FLOAT),
};

#ifdef DOXYGEN
#define mexFunction mexLoadAgilentTraces
#endif

/* [traces dirs] = mexLoadAgilentTraces(filename) */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  char * usage = 
    "usage: [traces mainh blockh] = mexLoadAgilentTraces(filename)\n";
  char measfile[512];
  vfile_info_t info;
  mwSize outdims[3];
  float *outdataR;
  float *outdataI;

  varian_printf = mexCommonsPrintfCallback;

  if (nrhs != 1 || nlhs != 3 || !mxIsChar(prhs[0]))
    mexErrMsgTxt(usage);

  if (mxGetString(prhs[0], measfile, sizeof(measfile)))
    mexErrMsgTxt(usage);

  if (varian_load_info(measfile, &info)) {
    mexWarnMsgTxt(measfile);
    mexErrMsgTxt("cant load fid main header\n");
  }

  if (info.main_header.np * info.main_header.ntraces * info.main_header.nblocks == 0) {
    mexErrMsgTxt("zero-sized data\n");
  }

  /* Construct outdata */
  outdims[0] = info.main_header.np;
  outdims[1] = info.main_header.ntraces;
  outdims[2] = info.main_header.nblocks;

  plhs[0] = mxCreateNumericArray(3, outdims, mxSINGLE_CLASS, mxCOMPLEX);
  outdataR = mxGetData(plhs[0]);
  outdataI = mxGetImagData(plhs[0]);

  /* load scandata */
  if (varian_load_fid(measfile,
                      info.main_header.ntraces * info.main_header.nblocks,
                      info.main_header.np,
                      outdataR, outdataI, 1, 0))
    mexErrMsgTxt("error loading agilent traces\n");

  plhs[1] = mexStructToArray(NARRAY(mheader_fields), mheader_fields, 1, &info.main_header);
  plhs[2] = mexStructToArray(NARRAY(bheader_fields), bheader_fields, info.main_header.nblocks, info.block_headers);

  varian_free_info(&info);
}
