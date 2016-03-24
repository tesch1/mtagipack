/*
 * Copyright 2010-2014 Regents of the University of Minnesota.
 * This file subject to terms of Creative Commons CC BY-NC 4.0, see LICENSE.html
 *
 * author(s): Michael Tesch (tesch1@gmail.com),
 */
/*! \file
 * \brief	Save agilent traces from an array.
 *
 * \section     Building
 *
 * Use mtmatpack_setup.m or mtagipack_setup.m
 *
 * \section     Example
\verbatim
>> mexSaveAgilentTraces
Error using mexSaveAgilentTraces
usage: mexSaveAgilentTraces(filename [, mainh, blockh])

>> 
>> [t, mainh, blockh] = mexLoadAgilentTraces('/home/tesch/data/data.fid');
>> mexSaveAgilentTraces('/home/tesch/data/data2.fid', t);
>> mexSaveAgilentTraces('/home/tesch/data/data3.fid', t, mainh, blockh);
>> 
\endverbatim

 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "varian.h"
#include "premexh.h"
#include "mex.h"

#include "mexCommons.h" // for mexStructArray()

/*
 * these are identical to what's in mexLoadAgilentTraces.c ...
 */
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
#define mexFunction mexSaveAgilentTraces
#endif

/* mexSaveAgilentTraces(filename) */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  char * usage = "usage: mexSaveAgilentTraces(filename, fids [, mainh, blockh])\n";
  char measfile[512];
  const mwSize * dims;
  float *outdataR;
  float *outdataI;

  varian_printf = mexCommonsPrintfCallback;

  if ((nrhs != 2 && nrhs != 4)
      || nlhs != 0 || !mxIsChar(prhs[0]) || !mxIsComplex(prhs[1]) 
      || !mxIsSingle(prhs[1]))
    mexErrMsgTxt(usage);

  if (mxGetString(prhs[0], measfile, sizeof(measfile)))
    mexErrMsgTxt(usage);

  dims = mxGetDimensions(prhs[1]);

  mwSize nblocks;
  if (3 == mxGetNumberOfDimensions(prhs[1]))
    nblocks = dims[2];
  else if (2 == mxGetNumberOfDimensions(prhs[1]))
    nblocks = 1;
  else
    mexErrMsgTxt(usage);

  outdataR = mxGetData(prhs[1]);
  outdataI = mxGetImagData(prhs[1]);

  mwSize nent = dims[0] * dims[1] * nblocks;
  float2 * traces = malloc(sizeof(float2) * nent);

  mwSize ii;
  for (ii = 0; ii < nent; ii++) {
    traces[ii].x = outdataR[ii];
    traces[ii].y = outdataI[ii];
  }

  /* save scandata */
  if (nrhs == 2)
    varian_save_fid(measfile, dims[1], dims[0], nblocks, traces);
  else {
    mwSize count;
    vmain_header_t * mainh = mexArrayToStruct(NARRAY(mheader_fields), mheader_fields, prhs[2], &count);
    vblock_header_t * blockh = mexArrayToStruct(NARRAY(bheader_fields), bheader_fields, prhs[3], &count);
    vfile_info_t vinfo;
    memcpy(&vinfo, mainh, sizeof(vmain_header_t));
    vinfo.main_header.np *= 2;
    vinfo.block_headers = blockh;
    varian_save_fids(measfile, &vinfo, traces);
    free(mainh);
    free(blockh);
  }
}
