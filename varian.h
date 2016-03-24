/*
 * Copyright 2010-2014 Regents of the University of Minnesota.
 * This file subject to terms of Creative Commons CC BY-NC 4.0, see LICENSE.html
 *
 * author(s): Michael Tesch (tesch1@gmail.com)
 *
 */
/*! \file
 * \brief	Interface to agilent/(once known as varian) io code.
 */
#ifndef VARIAN_H
#define VARIAN_H

#include <stdint.h>
#include "matvec.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * snooped the file format from
 * ftp://davinci.chem.umanitoba.ca/pub/marat/SpinWorks/source_library/Varian_header.cs
 */

/** varian fid file main header */
#pragma pack(push, 2)
typedef struct vfile_header_ {
  /*                    // byte index (offset) */
  int32_t nblocks;      /* 0 number of blocks */
  int32_t ntraces;      /* 4 traces per block */
  int32_t np;           /* 8 scalar samples per trace->complex samples per trace */
  int32_t ebytes;       /* 12 */
  int32_t tbytes;       /* 16 */
  int32_t bbytes;       /* 20 */
  uint16_t status;      /* 26 */
  uint16_t vers_id;     /* 24 */
  int32_t nbheaders;    /* 28 */
} vmain_header_t;

/** varian fid file block header */
typedef struct vblock_header_ {
  /*                    // byte index (offset) */
  uint16_t status;      /* 2 */
  uint16_t scale;       /* 0 */
  uint16_t mode;        /* 6 */
  uint16_t index;       /* 4 */
  int32_t ctcount;      /* 8 */
  float lpval;          /* 12 */
  float rpval;          /* 16 */
  float lvl;            /* 20 */
  float tlt;            /* 24 */
} vblock_header_t;
#pragma pack(pop)

typedef struct vfile_info_ {
  vmain_header_t main_header;
  vblock_header_t * block_headers;
  FILE * fp;
} vfile_info_t;

int varian_load_info(const char *dirname, vfile_info_t * vfile);
void varian_free_info(vfile_info_t * vfile);

int varian_load_mh(const char *dirname, vmain_header_t *mainh);
int varian_fid_ntraces(const char *dirname, size_t * ntracesP);
int varian_fid_np(const char *dirname, size_t * npP);
int varian_fid_nblocks(const char *dirname, size_t * nblocksP);
int varian_load_fid(const char *dirname, size_t ntraces, size_t np, float *re, float *im,
                    int nchan, int chan);
int varian_save_fid(const char *dirname, size_t ntraces, size_t np, int nblocks, const float2 *traces);
int varian_save_fids(const char *dirname, const vfile_info_t * vfile, const float2 *traces);
int varian_load_rf(const char *filename, size_t * npoints, vec2_t ** rf);

extern int (*varian_printf)(const char *format, ...);

#ifdef __cplusplus
}
#endif
#endif
