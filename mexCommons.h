/*
 * Copyright 2013 Regents of the University of Minnesota.
 * This file subject to terms of Creative Commons CC BY-NC 4.0, see LICENSE.html
 *
 * author(s): Michael Tesch (tesch1@gmail.com),
 */
#ifndef MEX_COMMONS_H
#define MEX_COMMONS_H
#include "mex.h"

#ifdef __cplusplus
extern "C" {
#endif

enum _ftype {TY_U32, TY_U16, TY_S16, TY_FLOAT, TY_S32, TY_DOUBLE};
size_t sizeof_ftype(enum _ftype ftype);

/*! descriptions of a field in a structure
 */
typedef struct {
  const char * name;    /*!<    field name */
  size_t offset;        /*!<    byte offset from start of struct */
  size_t size;          /*!<    field size in bytes */
  enum _ftype type;     /*!<    type of data stored in field */
} struct_fielddesc_t;

  mxArray * mexStructToArray(struct_fielddesc_t * sfields, size_t sizeofstruct, 
                             size_t count, void * data);
void * mexArrayToStruct(struct_fielddesc_t * sfields, size_t sizeofstruct, 
                        const mxArray * array, mwSize * count);

int mexCommonsPrintfCallback(const char *format, ...);

#define MC_FIELD_DEF(fname, ftype, structtype)                          \
  { #fname,                                                             \
      (size_t)&((structtype*)NULL)->fname,                              \
      sizeof(((structtype*)NULL)->fname),                               \
      ftype }

#ifndef MAX
#define MAX(a,b) ((a)>(b)?(a):(b))
#endif

#ifdef __cplusplus
}

#include <iostream>

/* class for redirecting cout to matlab's output stream */
class mexstream : public std::streambuf {
public:
  /* call at beginning of mexFunction() */
  static void install() {
    oldoutbuf(std::cout.rdbuf(&getmexout()));
  }
  static void restore() {
    std::cout.rdbuf(oldoutbuf());
  }
private:
  static mexstream & getmexout() {
    static mexstream mexout;
    return mexout;
  }
  static std::streambuf * oldoutbuf(std::streambuf * oldbuf = NULL) {
    static std::streambuf * oldoutbuf;
    std::streambuf * tmp;
    tmp = oldoutbuf;
    oldoutbuf = oldbuf;
    return tmp;
  }
  void drawnow() {
#if 0
#if 1
    /*mexEvalString("drawnow;"); /* let matlab flush the printf */
    /*mexEvalString("pause(0.0001);"); /* let matlab flush the printf */
    mexEvalString("drawnow;"); /* let matlab flush the printf */
#else
    mxArray * exception;
    exception = mexCallMATLABWithTrap(0, NULL, 0, NULL, "drawnow"); /* let matlab flush the printf */
    if(exception != NULL) {
      fprintf(stderr, "EXCEPTION\n");
      /* Throw the MException returned by mexCallMATLABWithTrap
       * after cleaning up any dynamically allocated resources */
      mexCallMATLAB(0, (mxArray **)NULL, 
                    1, &exception, "throw");
    }
#endif
#endif
  }
protected:
  virtual std::streamsize xsputn(const char *s, std::streamsize n) {
    mexPrintf("%.*s",n,s);
    drawnow();
    return n;
  }
  virtual int overflow(int c = EOF) {
    if (c != EOF) {
      mexPrintf("%.1s",&c);
      drawnow();
    }
    return 1;
  }
};

#endif
#endif /* MEX_COMMON_H */
