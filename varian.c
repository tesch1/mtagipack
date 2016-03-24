/*
 * Copyright 2010-2014 Regents of the University of Minnesota.
 * This file subject to terms of Creative Commons CC BY-NC 4.0, see LICENSE.html
 *
 * author(s): Michael Tesch (tesch1@gmail.com)
 *
 */
/*! \file
 * \brief Functions for reading varian .fid data blocks
 *
 * The test program in this file calculates some statistics on the sample data.
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <assert.h>
#include <arpa/inet.h>
#include <sys/stat.h>
#include <limits.h>
#include <ctype.h>
#include <errno.h>

#include "cprocpar.h"
#include "debug.h"
#include "varian.h"

/* header and block status bits */
#define V_STATUS_DATA          (0x0001)         //!< Data valid
#define V_STATUS_SPEC          (0x0002)         //!< Data is spectrum
#define V_STATUS_32BIT         (0x0004)         //!< Data is 32bit, else 16bit (valid if FLOAT==0)
#define V_STATUS_FLOAT         (0x0008)         //!< Data is 'float' not int
#define V_STATUS_COMPLEX       (0x0010)         //!< Data is 2-channel complex
#define V_STATUS_HCOMPLEX      (0x0020)         //!< Data is hypercomplex (??)

/* header only status bits */
#define VH_STATUS_DDR           (0x0080)        //!< Data acquired on DDR
#define VH_STATUS_2ND_FT        (0x0100)        //!< Second FT (huh)
#define VH_STATUS_TRANPOSED     (0x0200)        //!< Data transposed
#define VH_STATUS_3D            (0x0400)        //!< 3D data (huh)
#define VH_STATUS_NP            (0x0800)        //!< NP dim is active
#define VH_STATUS_NF            (0x1000)        //!< NF dim is active
#define VH_STATUS_NI            (0x2000)        //!< NI dim is active
#define VH_STATUS_NI2           (0x4000)        //!< NI2 dim is active

#define V_STATUS_IS_VALID(s)    (((unsigned int)(s) & V_STATUS_DATA) != 0)
#define V_STATUS_IS_INT16(s)    (((unsigned int)(s) & V_STATUS_FLOAT) == 0 && \
                                 ((unsigned int)(s) & V_STATUS_32BIT) == 0)
#define V_STATUS_IS_INT32(s)    (((unsigned int)(s) & V_STATUS_FLOAT) == 0 && \
                                 ((unsigned int)(s) & V_STATUS_32BIT) != 0)
#define V_STATUS_IS_FLOAT(s)    (((unsigned int)(s) & V_STATUS_FLOAT) != 0)
/* this complex flag is apparently not always valid */
#define xV_STATUS_IS_COMPLEX(s) (((unsigned int)(s) & V_STATUS_COMPLEX) == 1)
#define xV_STATUS_IS_HYPERCOMPLEX(s) (((unsigned int)(s) & V_STATUS_HCOMPLEX) == 1)
/* header only tests */
#define V_STATUS_IS_DDR(s)      (((unsigned int)(s) & VH_STATUS_DDR) != 0)      //!< ddr
#define V_STATUS_IS_TRANSPOSED(s) (((unsigned int)(s) & VH_STATUS_TRANPOSED) != 0) //!< transposed
#define V_STATUS_IS_3D(s)       (((unsigned int)(s) & VH_STATUS_3D) != 0)       //!< 3d
#define V_STATUS_IS_NP(s)       (((unsigned int)(s) & VH_STATUS_NP) != 0)       //!< np dim
#define V_STATUS_IS_NF(s)       (((unsigned int)(s) & VH_STATUS_NF) != 0)       //!< nf
#define V_STATUS_IS_NI(s)       (((unsigned int)(s) & VH_STATUS_NI) != 0)       //!< ni
#define V_STATUS_IS_NI2(s)      (((unsigned int)(s) & VH_STATUS_NI2) != 0)      //!< ni2

/*! \brief      module-local variable to remember if in verbose mode
 * 
 */
static int m_verbose = 0;
static int stderr_printf_callback(const char *format, ...);

//! change this to get error messages from varian.c functions
int (*varian_printf)(const char *format, ...) = stderr_printf_callback;

/* internal helper function */
static int stderr_printf_callback(const char *format, ...)
{
  va_list args;
  int ret;
  va_start(args, format);
  ret = vfprintf(stderr, format, args);
  va_end (args);
  return ret;
}

/* internal helper function */
static inline void swap4(void *d, size_t size)
{
  int32_t *data = d;
  size_t i;
  size /= sizeof(*data);
  for (i = 0; i < size; i++)
    data[i] = ntohl(data[i]);
}

/* internal helper function */
static inline void swap2(void *d, size_t size)
{
  int16_t *data = d;
  size_t i;
  size /= sizeof(*data);
  for (i = 0; i < size; i++)
    data[i] = ntohs(data[i]);
}

/* internal helper function */
static FILE * varian_fopen_fid(const char *dirname, const char *mode)
{
  char fidfilename[PATH_MAX];
  FILE *fp;

  /* check for dirname/fid */
  snprintf(fidfilename, sizeof(fidfilename), "%s/fid", dirname);
  fp = fopen(fidfilename, mode);
  if (fp)
    return fp;

  /* check for dirname/acqfil/fid */
  snprintf(fidfilename, sizeof(fidfilename), "%s/acqfil/fid", dirname);
  fp = fopen(fidfilename, mode);
  if (fp)
    return fp;

  /* give up */
  return NULL;
}

void print_mainh(vmain_header_t * mainh)
{
  varian_printf("mainh: nb:%d nt:%d np:%d status:0x%x vers_id:%d nbheaders:%d bbytes:%d\n",
         mainh->nblocks, mainh->ntraces, mainh->np, mainh->status, mainh->vers_id,
         mainh->nbheaders, mainh->bbytes);

  varian_printf("status flags: ");
  if (V_STATUS_IS_VALID(mainh->status)) varian_printf("valid,");
  if (V_STATUS_IS_INT16(mainh->status)) varian_printf("int16,");
  if (V_STATUS_IS_INT32(mainh->status)) varian_printf("int32,");
  if (xV_STATUS_IS_COMPLEX(mainh->status)) varian_printf("complex,");
  if (xV_STATUS_IS_HYPERCOMPLEX(mainh->status)) varian_printf("hypercomplex,");
  if (V_STATUS_IS_DDR(mainh->status)) varian_printf("ddr,");
  if (V_STATUS_IS_TRANSPOSED(mainh->status)) varian_printf("transposed,");
  if (V_STATUS_IS_3D(mainh->status)) varian_printf("3D,");
  if (V_STATUS_IS_NP(mainh->status)) varian_printf("np,");
  if (V_STATUS_IS_NF(mainh->status)) varian_printf("nf,");
  if (V_STATUS_IS_NI(mainh->status)) varian_printf("ni,");
  if (V_STATUS_IS_NI2(mainh->status)) varian_printf("ni2");
  varian_printf("\n");
}

/*
 * helper function
 */
int varian_load_mh(const char *dirname, vmain_header_t *mainh)
{
  FILE *fp;
  struct stat statbuf;

  fp = varian_fopen_fid(dirname, "r");
  if (!fp) {
    varian_printf("can't open fid file in %s\n", dirname);
    return -1;
  }

  /* read file header */
  if (1 != fread(mainh, sizeof(*mainh), 1, fp)) {
    fclose(fp);
    IF_FAIL(1, -1, "fread");
  }

  if (fstat(fileno(fp), &statbuf)) {
    varian_printf("unable to fstat fid in %s\n", dirname);
    fclose(fp);
    return -1;
  }
  fclose(fp);

  /* correct endianness */
  swap4(mainh, sizeof(*mainh));

  /* we consider np as count of complex points */
  mainh->np /= 2;

  /* verify file size against header info */
  if (statbuf.st_size != sizeof(vmain_header_t) + (size_t)mainh->nblocks * (size_t)mainh->bbytes) {
    varian_printf("fid file has inconsistent size %s (%zd != %zd)\n", dirname,
                  (long)statbuf.st_size,
                  sizeof(vmain_header_t) + (size_t)mainh->nblocks * (size_t)mainh->bbytes);
    print_mainh(mainh);
    return -1;
  }

  return 0;
}

int varian_load_info(const char *dirname, vfile_info_t * vfile)
{
  int ret;

  /* get main header */
  ret = varian_load_mh(dirname, &vfile->main_header);
  if (ret) {
    varian_printf("varian_load_info: unable to read master header\n");
    return ret;
  }

  /* open the file */
  vfile->fp = varian_fopen_fid(dirname, "r");
  if (!vfile->fp) {
    varian_printf("varian_load_info: unable to reopen fid file %s", dirname);
    return -1;
  }

  /* allocate memory for block headers */
  vfile->block_headers = calloc(vfile->main_header.nblocks, sizeof(vblock_header_t));
  if (!vfile->block_headers) {
    varian_printf("varian_load_info: calloc failed\n");
    varian_free_info(vfile);
    return -1;
  }

  /* read in block headers */
  off_t ii;
  for (ii = 0; ii < vfile->main_header.nblocks; ii++) {
    off_t boff = sizeof(vmain_header_t) + ii * vfile->main_header.bbytes;

    /* header offset */
    if (fseeko(vfile->fp, boff, SEEK_SET)) {
      print_mainh(&vfile->main_header);
      varian_printf("varian_load_info: fseek (%zd) failed\n", boff);
      varian_free_info(vfile);
      return -1;
    }

    /* read the header */
    if (1 != fread(&vfile->block_headers[ii], sizeof(vblock_header_t), 1, vfile->fp)) {
      varian_printf("varian_load_info: fread failed\n");
      varian_free_info(vfile);
      return -1;
    }

    /* correct block header endianness */
    swap4(&vfile->block_headers[ii], sizeof(vblock_header_t));
    
    /* check sanity? */
  }

  return 0;
}

/** \brief      Cleanup a vfile_info_t struct
 *
 * \param[in]	vfile   vfile_info_t struct to close down.
 */
void varian_free_info(vfile_info_t * vfile)
{
  if (vfile->block_headers) free(vfile->block_headers);
  if (vfile->fp) fclose(vfile->fp);
  vfile->block_headers = NULL;
  vfile->fp = NULL;
}

/** \brief      Get number of traces in a .fid directory's fid.
 *
 * \param[in]	dirname	name of varian format .fid directory
 */
int varian_fid_ntraces(const char *dirname, size_t * ntracesP)
{
  vmain_header_t mainh;
  IF_FAIL(varian_load_mh(dirname, &mainh), -1, "varian_fid_ntraces");
  *ntracesP = mainh.nblocks * mainh.ntraces;
  return 0;
}

/** \brief      Get number of complex data points per trace in a .fid directory's fid.
 *
 * \param[in]	dirname	name of varian format .fid directory
 */
int varian_fid_np(const char *dirname, size_t * npP)
{
  vmain_header_t mainh;
  IF_FAIL(varian_load_mh(dirname, &mainh), -1, "varian_fid_np");
  *npP = mainh.np;
  return 0;
}

/** \brief      Get number of data blocks in .fid/fid
 *
 * \param[in]	dirname	name of varian format .fid directory
 */
int varian_fid_nblocks(const char *dirname, size_t * nblocksP)
{
  vmain_header_t mainh;
  IF_FAIL(varian_load_mh(dirname, &mainh), -1, "varian_fid_nblocks");
  *nblocksP = mainh.nblocks;
  return 0;
}

/** \brief      Read a complex fid in from a varian (/agilent) .fid/ directory.
 *
 * \param[in]		dirname	path to the *.fid/ directory holding the /fid file
 * \param[in]		ntraces	number of traces to be expected in fid
 * \param[in]		np	number of complex points per trace
 * \param[in,out]	re	pointer to array to receive the real part of the data
 * \param[in,out]	im	pointer to array to receive the imag part of the data
 * \param[in]		nchan	number of channels in this dataset
 * \param[in]		chan	which channel to read
 *
 * re and im are both real arrays of dimension [ntraces * np].
 * \return      data pointer in samples, and the number of complex samples read. -1 for failure.
 */
int varian_load_fid(const char *dirname,
                    size_t ntraces, size_t np,
                    float *re, float *im,
                    int nchan, int chan)
{
  FILE *fp;
  vmain_header_t mainh;
  int b;
  size_t current;
  int32_t *blockdata;

  fp = varian_fopen_fid(dirname, "r");
  if (!fp) {
    varian_printf("varian_load_fid failed (%s)\n", dirname);
    return -1;
  }

  /* read file header to move file pointer */
  if (1 != fread(&mainh, sizeof(mainh), 1, fp)) {
    fclose(fp);
    varian_printf("varian_load_fid fread failed (%s)\n", dirname);
    return -1;
  }

  /* really read & decode header into local structure (hackish) */
  if (varian_load_mh(dirname, &mainh)) {
    varian_printf("varian_load_mh failed");
    return -1;
  }

  if (m_verbose) {
    print_mainh(&mainh);
  }

  /* some sanity checks */
  assert(ntraces == (size_t)(mainh.nblocks * mainh.ntraces) / nchan);
  assert(nchan > 0);
  assert(nchan <= 64); // arbitrary, but good sanity check for now
  assert(chan < nchan);
  assert(chan >= 0);
  assert(np == (size_t)mainh.np);

  // check block size
  {
    size_t samplebytes = V_STATUS_IS_INT16(mainh.status) ? sizeof(short) : sizeof(float);
    if ((size_t)mainh.bbytes - sizeof(vblock_header_t)
        != (mainh.np * 2 * mainh.ntraces * samplebytes)) {
      print_mainh(&mainh);
      assert(0);
    }
  }

  /* check block count */

  //assert((size_t)mainh.bbytes <= (mainh.np * 4 * mainh.ntraces * samplebytes));
  if (!V_STATUS_IS_VALID(mainh.status)) {
    varian_printf("no valid fid data? %x", mainh.status);
    return -1;
  }
  //IF_DIE(!V_STATUS_IS_COMPLEX(mainh.status), "fid data not complex? %x", mainh.status);

  /* check sample space */
  IF_FAIL(!re || !im, -1, "bad pointers");

  current = 0;

  /* allocate temporary memory for loading fid */
  blockdata = malloc(mainh.bbytes);
  IF_FAIL(!blockdata, -1, "malloc of %d bytes", mainh.bbytes);

  /* read block headers */
  for (b = 0; b < mainh.nblocks; b++) {
    vblock_header_t blockh;
    size_t blockbytes = mainh.bbytes - sizeof(blockh);
    int i;

    IF_FAIL(1 != fread(&blockh, sizeof(blockh), 1, fp), -1, "fread");

    /* correct endianness */
    swap4(&blockh, sizeof(blockh));

    if (m_verbose) {
      varian_printf("blockh: idx:%d status:%x scale:%d mode:%d ctcount:%d\n"
                    "               lpval:%f rpval:%f lvl:%f tlt:%f\n",
                    blockh.index, blockh.status, blockh.scale, blockh.mode, blockh.ctcount,
                    (double)blockh.lpval, (double)blockh.rpval, (double)blockh.lvl, (double)blockh.tlt);
    }

    /* check assumptions */
    IF_FAIL(!V_STATUS_IS_VALID(blockh.status), -1, 
            "no valid fid data in block %d? 0x%x", b, blockh.status);
    //IF_DIE(!V_STATUS_IS_COMPLEX(blockh.status), "fid data not complex? %x", mainh.status);

    /* read sample data */
    if (1 != fread(blockdata, blockbytes, 1, fp)) {
      perror("fread");
      varian_printf("fread b%d %ld", b, (long)blockbytes);
      return -1;
    }

    /*
     * decode each trace
     */
    int t;
    if (V_STATUS_IS_INT32(blockh.status)) {
      // 32 bit ints -- old-school inova
      //printf("32bit ints (idx %d)\n", blockh.index);
      int32_t *tracedata = blockdata + chan;
      int tracebytes = 4 * 2 * mainh.np;
      varian_printf("32bit ints (idx %d) untested!\n", blockh.index);
      for (t = 0; t < mainh.ntraces; t += nchan) {
        /* correct endianness */
        swap4(tracedata, tracebytes);
        for (i = 0; i < 2 * mainh.np; i += 2, current++) {
          re[current] = tracedata[i];
          im[current] = tracedata[i+1];
        }
        tracedata += 2 * mainh.np * nchan;
      }
    } else if (V_STATUS_IS_FLOAT(blockh.status)) {
      // 32 bit float, untested!!! TO DO
      //
      //varian_printf("32bit floats (idx %d) x %d\n", blockh.index, mainh.np * mainh.ntraces);
      float *tracedata = (float *)blockdata;
      int tracebytes = sizeof(float) * 2 * mainh.np;

      /* Curt says that the DDR data (which is saved as floats)
       * for multi-channel is every-other block, so see if we're
       * on the right block, otherwise continue
       */
      if ((b % nchan) != chan)
        continue;
      for (t = 0; t < mainh.ntraces; t++) {
        /* correct endianness */
        swap4(tracedata, tracebytes);
        for (i = 0; i < 2 * mainh.np; i += 2, current++) {
          re[current] = tracedata[i];
          im[current] = tracedata[i+1];
        }
        tracedata += 2 * mainh.np;
      }
    } else if (V_STATUS_IS_INT16(blockh.status)) {
      // 16 bit to save space (use STATUS!!!) TO DO
      // multi-channel data is 
      int16_t *tracedata = (int16_t *)blockdata;
      int tracebytes = sizeof(int16_t) * 2 * mainh.np;

      if ((b % nchan) != chan) // i think only on DDR?? should check ddr-status bit here
        continue;

      //varian_printf("16bit ints (idx %d) (untested!)\n", blockh.index);

      for (t = 0; t < mainh.ntraces; t++) {
        /* correct endianness */
        swap2(tracedata, tracebytes);
        for (i = 0; i < 2 * mainh.np; i += 2, current++) {
          re[current] = tracedata[i];
          im[current] = tracedata[i+1];
        }
        tracedata += 2 * mainh.np;
      }
    } else {
      IF_FAIL(1, -1, "unknown fid file block status type 0x%x", blockh.status);
    }
  }
  free(blockdata);
  fclose(fp);
  assert((double)current ==
         ((double)mainh.nblocks * (double)mainh.ntraces * (double)mainh.np) / (double)nchan);
  varian_printf("read blocks:%d traces:%d samples:%d channel:%d/%d\n", mainh.nblocks,
                mainh.nblocks * mainh.ntraces,
                mainh.nblocks * mainh.ntraces * mainh.np, chan, nchan);

  return 0;
}

/** \brief      Save a complex fid in from a varian (/agilent) .fid/ directory.
 *
 * \param[in]		dirname	path to the *.fid/ directory holding the /fid file
 * \param[in]		ntraces	number of traces to be in fid
 * \param[in]		np	number of complex points per trace
 * \param[in,out]	re	pointer to array of real part of the data
 * \param[in,out]	im	pointer to array of imag part of the data
 *
 * re and im are both real arrays of dimension [ntraces * np].
 */
int varian_save_fids(const char *dirname, const vfile_info_t * vfile, const float2 *traces)
{
  FILE *fp;
  size_t i;
  int b;
  size_t blockbytes;
  float * blockdata;
  size_t np = vfile->main_header.np / 2;
  size_t ntraces = vfile->main_header.ntraces;
  size_t nblocks = vfile->main_header.nblocks;

  /* check if dirname exists, if so, fail */
  //struct stat statbuf;

  //if (!stat(dirname, &statbuf)) {
  //  printf("didn't save fid to '%s/fid' because %s already exists\n", dirname, dirname);
  //  return;
  //}

  /* create dirname as directory */
  if (mkdir(dirname, 0700) && errno != EEXIST) {
    varian_printf("unable to mkdir(%s), not saving fid\n", dirname);
    return -1;
  }

  /* open dirname/fid for output */
  fp = varian_fopen_fid(dirname, "w");
  if (!fp) {
    varian_printf("can't open fid file in %s", dirname);
    return -1;
  }


  /* write the main header */
  {
    vmain_header_t mainh;

    memcpy(&mainh, &vfile->main_header, sizeof(vmain_header_t));

    swap4(&mainh, sizeof(vmain_header_t));

    if (1 != fwrite(&mainh, sizeof(mainh), 1, fp)) {
      varian_printf("varian_save_fids() fwrite failed");
      return -1;
    }
  }

  blockbytes = sizeof(float) * 2 * np * ntraces;
  blockdata = malloc(blockbytes);
  if (!blockdata) {
    varian_printf("varian_save_fids malloc failed %ld", blockbytes);
    return -1;
  }

  for (b = 0; b < nblocks; b++) {
    const float2 * block;

#if 0
    varian_printf("blockh: idx:%d status:%x scale:%d mode:%d ctcount:%d\n"
                  "               lpval:%f rpval:%f lvl:%f tlt:%f\n",
                  blockh.index, blockh.status, blockh.scale, blockh.mode, blockh.ctcount,
                  (double)blockh.lpval, (double)blockh.rpval, (double)blockh.lvl, (double)blockh.tlt);
#endif

    /* correct endianness */
    {
      vblock_header_t blockh;

      memcpy(&blockh, &vfile->block_headers[b], sizeof(vblock_header_t));

      // 32 bit float, untested!!! TO DO
      varian_printf("32bit floats (idx %d)\n", blockh.index);

      swap4(&blockh, sizeof(blockh));

      if (1 != fwrite(&blockh, sizeof(blockh), 1, fp)) {
        varian_printf("varian_save_fids() : write block header failed");
        return -1;
      }
    }

    /* write block */
    block = &traces[b * ntraces * np];

    for (i = 0; i < np * ntraces; i++) {
      blockdata[2*i    ] = block[i].x;
      blockdata[2*i + 1] = block[i].y;
    }
    swap4(blockdata, blockbytes);

    if (1 != fwrite(blockdata, blockbytes, 1, fp)) {
      varian_printf("varian_save_fids() fwrite failed2");
      return -1;
    }
  }

  free(blockdata);
  fclose(fp);

  varian_printf("wrote fid (%d x %d) x %d blocks to %s\n", 
                (int)ntraces, (int)np, nblocks, dirname);

  return 0;
}

int varian_save_fid(const char *dirname, size_t ntraces, size_t np, int nblocks, const float2 * traces)
{
  vfile_info_t vfile;

  bzero(&vfile, sizeof(vfile));

  vfile.block_headers = calloc(nblocks, sizeof(vblock_header_t));
  if (!vfile.block_headers) {
    varian_free_info(&vfile);
    return -1;
  }

  /* build file header to move file pointer */
  vfile.main_header.nblocks = nblocks;
  vfile.main_header.ntraces = ntraces;
  vfile.main_header.np = np * 2; // complex data
  vfile.main_header.ebytes = sizeof(float);
  vfile.main_header.tbytes = 0;
  vfile.main_header.bbytes = sizeof(float) * 2 * np * ntraces + sizeof(vblock_header_t);
  vfile.main_header.status = V_STATUS_DATA | V_STATUS_FLOAT /*| V_STATUS_0x40?*/;
  vfile.main_header.vers_id = 0;
  vfile.main_header.nbheaders = 1;

#if 0
  varian_printf("write mainh: nblocks:%d ntraces:%d np:%d status:0x%x vers_id:%d nbheaders:%d bbytes:%d\n",
                vfile.main_header.nblocks, vfile.main_header.ntraces, vfile.main_header.np,
                vfile.main_header.status, vfile.main_header.vers_id,
                vfile.main_header.nbheaders, vfile.main_header.bbytes);
#endif

  /* fill-in block headers */
  int bb;
  for (bb = 0; bb < nblocks; bb++) {
    vfile.block_headers[bb].index = 1 + bb;
    vfile.block_headers[bb].status = V_STATUS_DATA | V_STATUS_FLOAT /*| V_STATUS_0x40?*/;
    vfile.block_headers[bb].scale = 1;
    vfile.block_headers[bb].mode = 0;
    vfile.block_headers[bb].ctcount = 1;
    vfile.block_headers[bb].lpval = 1.0;
    vfile.block_headers[bb].rpval = 1.0;
    vfile.block_headers[bb].lvl = 1.0;
    vfile.block_headers[bb].tlt = 0.0;
  }

  int ret = varian_save_fids(dirname, &vfile, traces);
  varian_free_info(&vfile);
  return ret;
}

#ifdef TEST

#include "array.h"

/** \brief	calculate & print some statistics on an array of reals
 *
 * \param[in]	n	array length
 * \param[in]	a	pointer to array
 */
void dist_stats(size_t n, real_t *a)
{
  size_t i, j, v, l;

  /* sort the array by value */
  array1r_sort(n, a);

  // count unique values
  v = 0;
  for (i = 1; i < n; i++)
    if (a[i - 1] != a[i])
      v++;
  varian_printf("unique: %d = %2.2f %%\n", (int)v, 100.0 * (double)v / (double)n);

  // todo: calculate some % of used bits within the bit range of min-max

  // find longest repeating series (mode), and its length
  v = 0;
  l = 0;
  for (i = 0; i < n - 1; i++) {
    j = i + 1;
    while (a[i] == a[j] && i < n && j < n)
      j++;
    if (j - i > v) {
      v = j - i;
      l = i;
    }
    i = j - 1;
  }
  varian_printf("longest sequence (mode): val:%12.12f (len: %d %2.2f %%)\n", 
                a[l], (int)v, 100.0 * (double) v / (double)n);
}

int varian_load_rf(const char *filename, size_t * npoints, vec2_t ** rf)
{
  // it's not, load from file
  FILE *fp;
  int isDEC = 0;
  int linenum = 0;
  char line[1024];
  size_t len;
  int cnt;
  vec2_t *c = NULL;
  double max = 0;

  fp = fopen(filename, "r");
  if (!fp)
    return -1;

  if (strstr(filename, ".DEC"))
    isDEC = 1;

  len = 0;
  while (fgets(line, sizeof(line), fp)) {
    double mag, phs, dur;
    int gate;
    linenum++;

    if (line[0] == '#')
      continue;

    // check if empty line, then ignore
    size_t j;
    for (j = 0; j < strlen(line); j++)
      if (!isspace(line[j]))
        break;
    if (j == strlen(line))
      continue;

    cnt = sscanf(line, "%lf %lf %lf %d", &phs, &mag, &dur, &gate);
    if (cnt < 3) {
      fprintf(stderr, "error reading '%s' on line %d\n", filename, linenum);
      fclose(fp);
      if (c) free(c);
      return -1;
    }

    if (isDEC) {
      double tmp_phs, tmp_mag, tmp_dur;
      tmp_phs = mag;
      tmp_mag = dur;
      tmp_dur = phs;
      gate = gate != 0;
      phs = tmp_phs;
      mag = tmp_mag;
      dur = (int)(tmp_dur / 90.0);
    }

    if (cnt == 2)
      dur = 1;

    // scale param -- sometimes there will be an entry of dur=0 or gate=0
    if (mag > max)
      max = mag;

    // gate is simple, for now
    if (cnt < 4)
      gate = 1;
    else
      if (!gate)
        mag = 0;

    // check bounds
    if (mag < 0.0) {
      printf("error in %s, line %d, mag out of bounds %f\n", filename, linenum, mag);
      mag = 0;
    }
    if (dur > (1<<16)) {
      printf("error in %s, line %d, duration too big: %f\n", filename, linenum, dur);
      fclose(fp);
      if (c) free(c);
      return -1;
    }

    // expand dur
    while (dur > 0.0) {
      dur -= 1.0;

      // realloc and insert
      c = (vec2_t *)realloc(c, (len + 1) * sizeof(vec2_t));
      IF_DIE(!c, "malloc fail");

      // append new value
      c[len].x = creal(mag * cexp(-I * deg_to_rad(phs)));
      c[len].y = cimag(mag * cexp(-I * deg_to_rad(phs)));

      len++;
    }
  }

  // rescale
  for (size_t i = 0; i < len; i++) {
    c[i].x *= 1. / max;
    c[i].y *= 1. / max;
  }

  if (!len) {
    fprintf(stderr, "pulse file empty: '%s'", filename);
    fclose(fp);
    if (c) free(c);
    return -1;
  }

  // done with this
  fclose(fp);

  // return the pulse
  *npoints = len;
  *rf = c;

  return 0;
}

#include "cargs.h"

static argspec_t opts[] = {
  {'t', "trace",        "trace index to print", ARG_INT},
  {'o', "outfile",      "filename to save traces to", ARG_STRING},
  {'v', "verbose",      "print info when reading headers", ARG_FLAG},
};

/** \brief      test the varian data loader and run some statistics on the trace data.
 *
 */
int main(int argc, char *argv[])
{
  real_t *re, *im;
  size_t ns, np, ntraces;
  char * trace_filename = NULL;
  int t = 0;
  int nch = 1;
  int ch;

  char *dirname = argv[1];

  int i;
  argspec_t *opt = NULL;
  for (i = 0; i < argc; i++) {
    opt = cgetopt(argc, argv, NARRAY(opts), opts);
    if (!opt) {
      break;
    }
    //varian_printf("arg %d: '%c'\n", i, opt->flag);
    switch (opt->flag) {
    case 'o': // output filename
      trace_filename = opt->val.s;
      break;
    case 't': // which trace to print
      t = opt->val.i;
      break;
    case 'v': // verbose output
      m_verbose = 1;
      break;
    }
  }

  /* try to open procpar if available */
  procpar_t * procpar;
  char parpath[PATH_MAX];
  snprintf(parpath, sizeof(parpath), "%s/procpar", dirname);
  procpar = pp_open(parpath, NULL, TREE_CURRENT);
  if (procpar) {
    const char * rcvrs;
    pp_get_string(procpar, "rcvrs", 0, &rcvrs);
    if (rcvrs) {
      nch = 0;
      while (*rcvrs)
        if (*rcvrs++ == 'y')
          nch++;
      varian_printf("found procpar, detected %d channels\n", nch);
    }
    pp_close(procpar);
  }

  /* load the fid data */
  varian_printf("loading fid\n");
  if (varian_fid_ntraces(dirname, &ntraces) || varian_fid_np(dirname, &np))
    return -1;
  assert(!(ntraces % nch));
  ntraces /= nch;
  ns = np * ntraces;

  for (ch = 0; ch < nch; ch++) {
    {
      float *RE, *IM;
      RE = calloc(ns, sizeof(float));
      IM = calloc(ns, sizeof(float));

      if (varian_load_fid(dirname, ntraces, np, RE, IM, nch, ch)) {
        varian_printf("unable to load fid");
        return -1;
      }
      varian_printf("ntraces:%d np:%d ns:%d\n", (int)ntraces, (int)np, (int)ns);

      re = array1r_clonef(ns, RE);
      im = array1r_clonef(ns, IM);
      free(RE);
      free(IM);
    }

    if (trace_filename) {
      complex_t * tmp = array1c_alloc(np);
      size_t idx = t * np;
      array_r_to_c(np, tmp, &re[idx], &im[idx]);
      array1c_txt(trace_filename, np, tmp);
      array1c_free(tmp);
    }

    /*
     * run stats
     */

    real_t min, max, avg, stddev;

    /* calculate some stats on the sample space */
    array1r_stats(ns, re, &min, &max, &avg, &stddev);
    varian_printf("real: min:%f max:%f avg:%f stddev:%f\n", min, max, avg, stddev);
    array1r_stats(ns, im, &min, &max, &avg, &stddev);
    varian_printf("imag: min:%f max:%f avg:%f stddev:%f\n", min, max, avg, stddev);

    /* try to find some kind of resolution measure for the samples */
    varian_printf("real:");
    dist_stats(ns, re);
    //array1r_txt("sorted_re.txt", ns, re);

    varian_printf("imag:");
    dist_stats(ns, im);
    //array1r_txt("sorted_im.txt", ns, im);

    /* cleanup */
    free(re);
    free(im);
  }
  return 0;
}

#endif
