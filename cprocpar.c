/*
 * Copyright 2010-2014 Regents of the University of Minnesota.
 * This file subject to terms of Creative Commons CC BY-NC 4.0, see LICENSE.html
 *
 * author(s): Michael Tesch (tesch1@gmail.com)
 *
 */
/*! \file
 * \brief utility functions for varian procpar formatted files
 */
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <ctype.h>
#include <errno.h>
#include <string.h>
#include <assert.h>
#include <regex.h>
#include <sys/stat.h>

#include "cprocpar.h"

#define NARRAY(a)       ((sizeof a) / sizeof (a)[0])
#define MAX_STR         (1024)          //!< many string statements limited by this

/* internal bookkeeping per-variable */
typedef struct ppvar_ {
  char name[128];
  enum ppvar_subtype subtype;
  enum ppvar_type type;
  enum ppvar_treemask tree;
  double min;
  double max;
  double step;
  int group;
  int display_group;
  int prot;
  int active;
  int nvals;
  double *realvalues;
  char **strvalues;
  int nenums;
  double *realenums;
  char **strenums;
  struct ppvar_ *next;
} ppvar_t;

/* internal bookkeeping per-procpar file */
struct procpar_ {
  ppvar_t *vars;
  ppvar_t *iter; /* for pp_get_next() */
};
//typedef struct procpar_ procpar_t;

/* bad global */
static int pp_line = 0;

//#if !defined(HAVE_STRNDUP) && !defined(strndup)
#if 0
static char * strndup(const char *str, size_t n)
{
  int len = strlen(str);
  len = n < len ? n : len;
  char *s = malloc(len + 1);
  assert(s);
  int i;
  for (i = 0; i < len; i++)
    s[i] = str[i];
  s[i] = 0;
  return s;
}
#endif

/*! \brief      Error-checked strtol (string to long int)
 *
 * Call strtol, and then check all the error conditions in one place.
 */
static long int strtol_x(const char *nptr, int base, int *err)
{
  char *endptr;
  errno = 0;
  long int i = strtol(nptr, &endptr, base);

  /* Check for various possible errors */
  if ((errno == ERANGE && (i == LONG_MAX || i == LONG_MIN))
      || (errno != 0 && i == 0)
      || (endptr && *endptr && !isspace(*endptr))) {
    if (errno && i == 0) perror("strtol");
    *err = ERANGE;
  }
  else if (endptr == nptr) {
    *err = ERANGE;
  }
  else
    *err = 0;
  return i;
}

/* internal function for reading array of n reals */
static double *pp_readreals(FILE *fp, int *n, enum ppvar_subtype subtype)
{
  double *vals;
  int i;

  // read the number of values in this variable
  if (1 != fscanf(fp, "%d", n)) {
    fprintf(stderr, "%s:%d fscanf failed\n", __FILE__, __LINE__);
    return NULL;
  }

  // allocate space for the vars
  vals = malloc(sizeof(double *) * *n);
  if (!vals) {
    fprintf(stderr, "%s:%d malloc failed\n", __FILE__, __LINE__);
    return NULL;
  }

  // read the variables
  for (i = 0; i < *n; i++) {
    double value;
    if (1 != fscanf(fp, "%lg", &value)) {
      fprintf(stderr, "error reading procpar line %d, invalid value",
              pp_line);
      return vals;
    }
    switch (subtype) {
    case UVAR_REAL:
    case UVAR_DELAY:
    case UVAR_FREQ:
      vals[i] = value;
      break;
    case UVAR_PULSE:
      vals[i] = value * 1e-6;
      break;
    case UVAR_INT:
      vals[i] = (double)(long long)value;
      break;
    case UVAR_STRING:
    case UVAR_FLAG:
    default:
      fprintf(stderr, "pp_readreals(): line %d, bad subtype %d\n",
              pp_line, subtype);
      vals[i] = value;
    }
  }

  return vals;
}

/* helper function for reading quoted strings */
static char * pp_readstring(FILE *fp)
{
  const int n = 1024; // max str len
  char *s;
  int i = 0;
  int ch, lastch;

  s = malloc(n);
  if (!s) {
    printf("%s:%d malloc failed\n", __FILE__, __LINE__);
    exit(-1);
  }

  // read until first "
  do {
    ch = fgetc(fp);
    if (EOF == ch) {
      printf("%s:%d fgetc failed @ line %d\n", __FILE__, __LINE__, pp_line);
      exit(-1);
    }
  } while (ch != '"');

  // read until next "
  do {
    lastch = ch;
    ch = fgetc(fp);
    if (EOF == ch || i >= n) {
      printf("%s:%d fgetc failed @ line %d\n", __FILE__, __LINE__, pp_line);
      exit(-1);
    }
    s[i] = (char)ch;
    i++;
  } while (!(ch == '"' && lastch != '\\'));

  // terminate string over the "
  s[i-1] = 0;

  return s;
}

/* internal function for reading array of n strings */
static char **pp_readstrs(FILE *fp, int *n)
{
  char **vals;
  int i;

  // read the number of values in this variable
  if (1 != fscanf(fp, "%d", n)) {
    fprintf(stderr, "%s:%d fscanf failed @~ line %d\n",
            __FILE__, __LINE__, pp_line);
    return NULL;
  }
  // allocate space for the vars
  vals = malloc(sizeof(*vals) * *n);
  if (!vals) {
    fprintf(stderr, "%s:%d malloc failed\n", __FILE__, __LINE__);
    return NULL;
  }

  // read the variables
  for (i = 0; i < *n; i++) {
    vals[i] = pp_readstring(fp);
    // check flag variables and enums
    //switch (subtype) {
    //case UVAR_STRING:
    //case UVAR_FLAG:
    //  break;
    //}
  }
  return vals;
}

/* internal function for making new procpar variables
 *
 * return NULL on error
 */
static ppvar_t * pp_newvar(const char *varname, 
                           enum ppvar_type type,
                           enum ppvar_subtype subtype,
                           enum ppvar_treemask tree)
{
  ppvar_t *v;

  if (strlen(varname) >= sizeof(v->name)) {
    fprintf(stderr, "pp_newvar(): variable name too long\n");
    return NULL;
  }
  v = calloc(1, sizeof(ppvar_t));
  if (!v) {
    fprintf(stderr, "out of memory, calloc() failed\n");
    return NULL;
  }
  strncpy(v->name, varname, sizeof(v->name));
  v->type = type;
  v->subtype = subtype;
  v->tree = tree;
  v->max = 1.e18;
  v->min = -1.e18;
  v->step = 0;
  v->nvals = 0;
  v->nenums = 0;
  v->active = 1;
  if (v->type == VAR_STRING) {
    v->strvalues = malloc(sizeof(char **));
    v->strenums = malloc(sizeof(char **));
  }
  else {
    v->realvalues = malloc(sizeof(double *));
    v->realenums = malloc(sizeof(double *));
  }
  v->next = NULL;
  return v;
}

/* 
 * the reverse of pp_newvar()
 */
static void pp_delvar(ppvar_t *v)
{
  assert(v);
  if (v->strvalues) free(v->strvalues);
  if (v->strenums) free(v->strenums);
  if (v->realvalues) free(v->realvalues);
  if (v->realenums) free(v->realenums);
  free(v);
}

/* internal function for adding new variable to a procpar set
 *
 */
static void pp_addvar(procpar_t *procpar, ppvar_t *var)
{
  var->next = procpar->vars;
  procpar->vars = var;
}

/* internal function for reading single procpar variables
 *
 * return NULL on errors
 */
static ppvar_t * pp_readvar(FILE *fp, enum ppvar_treemask tree)
{
  // read the info line
  int stype, type;
  int group, display_group, prot, active, i4;
  int res;
  char vname[140];
  double max, min, round;

  errno = 0;
  if (11 != (res = fscanf(fp, "%130s %d %d %lg %lg %lg %d %d %d %d %d",
                          vname, &stype, &type,
                          &max, &min, &round,
                          &group, &display_group, &prot, &active, &i4))) { // not sure these are right!!!
    // end of file condition
    if (res == EOF && !errno) {
      return NULL;
    }
    pp_line++;
    fprintf(stderr, "pp_readvar() fscanf failed, line %d, res=%d %p %s\n",
            pp_line, res, fp, vname);
    return NULL;
  }
  if (i4 != 64)
    fprintf(stderr, "pp_readvar() warning:%d int size(%d) != 64\n", pp_line, i4);
  if (active != 1 && active != 0)
    fprintf(stderr, "pp_readvar() warning:%d active == %d\n", pp_line, active);

  ppvar_t *v;
  v = pp_newvar(vname, type, stype, tree);
  if (!v)
    return NULL;
  v->max = max;
  v->min = min;
  v->step = round;
  v->group = group;
  v->display_group = display_group;
  v->prot = prot;
  v->active = active;

  switch (v->type) {
  case VAR_INT:
  case VAR_REAL:
    // read the variables
    free(v->realvalues);
    v->realvalues = pp_readreals(fp, &v->nvals, v->subtype);
    pp_line++;
    // read the enums
    free(v->realenums);
    v->realenums = pp_readreals(fp, &v->nenums, v->subtype);
    if (!v->realvalues || !v->realenums) {
      fprintf(stderr, "warning: procpar BAD VARIABLE '%s' ignoring rest of file %d\n",
              v->name, v->subtype);
      pp_delvar(v);
      return NULL;
    }
    pp_line++;
    break;
  case VAR_STRING:
    // read the variables
    free(v->strvalues);
    v->strvalues = pp_readstrs(fp, &v->nvals);
    if (!v->strvalues) {
      fprintf(stderr, "warning: procpar BAD VARIABLE '%s' ignoring rest of file %d\n",
              v->name, v->subtype);
      pp_delvar(v);
      return NULL;
    }
    pp_line += v->nvals;
    // read the enums
    free(v->strenums);
    v->strenums = pp_readstrs(fp, &v->nenums);
    if (!v->strenums) {
      fprintf(stderr, "warning: procpar BAD VARIABLE '%s' ignoring rest of file %d\n",
              v->name, v->subtype);
      pp_delvar(v);
      return NULL;
    }
    pp_line++;
    if (v->subtype != UVAR_STRING &&
        v->subtype != UVAR_FLAG) {
      fprintf(stderr, "warning: procpar BAD VARIABLE '%s' ignoring %d\n",
              v->name, v->subtype);
      v->active = 0;
    }
    break;
  case VAR_NOEXIST:
  default:
    fprintf(stderr, "pp_readvar() failed, line %d, vtype=%d %p vname=%s\n",
            pp_line, v->type, fp, v->name);
    pp_delvar(v);
    return NULL;
  }

  return v;
}

/*! \brief	create a new empty procpar. */
procpar_t * pp_new()
{
  procpar_t *pp = malloc(sizeof(procpar_t));
  if (!pp)
    return NULL;

  /* empty */
  pp->vars = NULL;

  /* initialize the iterator */
  pp->iter = NULL;

  return pp;
}

/*! \brief	open procpar, read it, and make it available.
 *
 * \param[in]   filename        file to read parameters from
 * \param[in]   old             if non-NULL, a procpar_t to add the variables from 'filename' to
 *                              existing vars in old are overwritten by values from 'filename'
 *
 * \return      0 for success, -1 for error
 */
procpar_t * pp_open(const char *filename, procpar_t *old, enum ppvar_treemask tree)
{
  ppvar_t *nvar;
  FILE *fp;
  procpar_t *pp;

  struct stat statbuf;
  if (stat(filename, &statbuf)) {
    fprintf(stderr, "unable to fstat %s\n", filename);
    return NULL;
  }
  if (S_ISREG(statbuf.st_mode)) {
  fp = fopen(filename, "r");
  if (!fp) {
    //fprintf(stderr, "unable to open %s\n", filename);
    return NULL;
  }
  } else {
    char procparfilename[PATH_MAX];
    snprintf(procparfilename, sizeof(procparfilename), "%s/procpar", filename);
    fp = fopen(procparfilename, "r");
    if (!fp) {
      //fprintf(stderr, "unable to open %s\n", filename);
      return NULL;
    }
  }

  pp = malloc(sizeof(procpar_t));
  if (!pp) {
    fclose(fp);
    return NULL;
  }
  pp->vars = NULL;

  //printf("opened '%s'  %p\n", filename, fp);
  pp_line = 0;

  /* read all the variables */
  while ((nvar = pp_readvar(fp, tree)))
    pp_addvar(pp, nvar);

  if (!feof(fp) || ferror(fp)) {
    fprintf(stderr, "pp_open: error reading %s around line %d\n", filename, pp_line);
  }

  fclose(fp);

  //printf("loaded '%s'  %d\n", filename, pp_size(pp));

  /* initialize the iterator */
  pp->iter = NULL;

  /* merge old into it */
  if (old) {
    pp = pp_merge(old, pp);
  }

  return pp;
}

/*! \brief     merge one procpar contents into another, close src
 *
 * \param[in]   dest    a parameter set to add the variables from 'filename' to
 *                      existing vars in dest are overwritten by values from src
 * \param[in]   src     a parameter set that is to be merged into dest
 *
 * \return      dest    returns dest
 */
procpar_t * pp_merge(procpar_t * dest, procpar_t * src)
{
  ppvar_t * v;

  //printf("merging  %d  %d\n", pp_size(dest), pp_size(src));

  /* delete duplicates in dest */
  for (v = src->vars; v && v->next; v = v->next)
    pp_delete(dest, v->name);

  //printf("merging  %d  %d\n", pp_size(dest), pp_size(src));

  /* splice src into dest */
  if (v) {
    pp_delete(dest, v->name);
    v->next = dest->vars;
    dest->vars = src->vars;
  }
  free(src);

  //printf("merged  %d\n", pp_size(dest));

  return dest;
}

/*! \brief     count how many variables in a parameter set */
int pp_size(procpar_t *procpar)
{
  int i = 0;
  ppvar_t *v;
  for (v = procpar->vars; v; v = v->next)
    i++;
  return i;
}

/* helper, frees the var, returns the next ptr */
static ppvar_t *pp_freevar(ppvar_t *v)
{
  ppvar_t *tmp;
  int i;
  if (v->type == VAR_STRING) {
    for (i = 0; i < v->nvals; i++)
      free(v->strvalues[i]);
    for (i = 0; i < v->nenums; i++)
      free(v->strenums[i]);
    free(v->strvalues);
    free(v->strenums);
  }
  else {
    free(v->realvalues);
    free(v->realenums);
  }
  tmp = v->next;
  free(v);
  return tmp;
}

/*! \brief	close procpar when no longer needed, deallocate memory */
void pp_close(procpar_t *procpar)
{
  ppvar_t *v;
  for (v = procpar->vars; v; )
    v = pp_freevar(v);
  free(procpar);
}

/* internal function to find a ppvar_t in the procpar_t matching 'varname' */
static ppvar_t * pp_findvar(const procpar_t *procpar, const char *varname)
{
  ppvar_t *v = procpar->vars;
  while (v && strcmp(v->name, varname))
    v = v->next;
  return v;
}

/*! \brief	delete a variable */
int pp_delete(procpar_t *procpar, const char * varname)
{
  ppvar_t *v;
  v = pp_findvar(procpar, varname);
  if (!v) {
    //fprintf(stderr, "unable to delete %s\n", varname);
    return -1;
  }

  //fprintf(stderr, " delete %s\n", varname);

  // point the prev at the next
  if (procpar->vars == v) {
    procpar->vars = v->next;
  } else {
    // find the previous in the chain
    ppvar_t * prev;
    for (prev = procpar->vars; prev && prev->next != v; prev = prev->next)
      ;
    assert(prev && prev->next == v); // bugs?
    prev->next = v->next;
  }

  // delete the var
  pp_freevar(v);

  return 0;
}

/*! \brief      Save params from specified tree(s) to a file. */
int pp_fwrite(FILE *file, procpar_t *procpar, enum ppvar_treemask tree)
{
  const char *varname;

  while ((varname = pp_get_next(procpar))) {
    ppvar_t * v;
    v = pp_findvar(procpar, varname);
    if (!(v->tree & tree))
      continue;
    fprintf(file, "%s %d %d %lg %lg %lg %d %d %d %d %d\n",
            varname,
            v->subtype, v->type,
            v->max, v->min, v->step,
            v->group, v->display_group, v->prot, v->active,
            (int)(8 * sizeof(long int)));
    fprintf(file, "%d ", v->nvals);
    int i;
    for (i = 0; i < v->nvals; i++) {
      switch (v->type) {
      case VAR_INT:
        fprintf(file, "%ld ", (long int)v->realvalues[i]);
        break;
      case VAR_REAL:
        fprintf(file, "%.12g ", v->realvalues[i]);
        break;
      case VAR_STRING:
        if (i) fprintf(file, "\n");
        fprintf(file, "\"%s\"", v->strvalues[i]);
        break;
      default:
        return -1;
      }
    }
    fprintf(file, "\n");
    fprintf(file, "%d ", v->nenums);
    for (i = 0; i < v->nenums; i++) {
      switch (v->type) {
      case VAR_INT:
        fprintf(file, "%ld ", (long int)v->realenums[i]);
        break;
      case VAR_REAL:
        fprintf(file, "%.12g ", v->realenums[i]);
        break;
      case VAR_STRING:
        fprintf(file, "\"%s\" ", v->strenums[i]);
        break;
      default:
        return -1;
      }
    }
    fprintf(file, "\n");
  }

/*
  if (11 != (res = fscanf(fp, "%130s %d %d %lg %lg %lg %d %d %d %d %d",
  vname, &stype, &type,
  &max, &min, &round,
  &i1, &i2, &i3, &i4, &i5))) {
*/
  return 0;
}

/*! \brief      Save procpar to a file. */
int pp_write(const char *filename, procpar_t *procpar, enum ppvar_treemask tree)
{
  FILE *fp;
  int ret;

  fp = fopen(filename, "w");
  if (!fp)
    return -1;

  ret = pp_fwrite(fp, procpar, tree);

  fclose(fp);
  return ret;
}

/*! \brief	print a variable */
void pp_print(procpar_t *procpar, const char *varname)
{
  int i;
  ppvar_t *v = pp_findvar(procpar, varname);
  if (!v) {
    printf("no such variable: %s\n", varname);
    return;
  }

  printf("%s %d : ", v->name, v->nvals);
  for (i = 0; i < v->nvals; i++)
    if (v->type == VAR_STRING)
      printf("\"%s\" ", v->strvalues[i]);
    else
      printf("%.8g ", v->realvalues[i]);
  printf("\n");
}

/*! \brief	iterates through the procpar's variable names,
 *              returning NULL at the turnover point.
 */
const char * pp_get_next(procpar_t *procpar)
{
  const char *varname;

  if (!procpar->iter) {
    procpar->iter = procpar->vars;
    varname = procpar->iter->name;
  } else {
    procpar->iter = procpar->iter->next;
    varname = procpar->iter ? procpar->iter->name : NULL;
  }

  return varname;
}

/*! \brief	inquire a particular variable's type
 */
enum ppvar_type pp_get_type(const procpar_t *procpar, const char *varname)
{
  ppvar_t *v = pp_findvar(procpar, varname);
  if (!v) return VAR_NOEXIST;
  return v->type;
}

/*! \brief	inquire a particular variable's subtype
 */
enum ppvar_subtype pp_get_subtype(const procpar_t *procpar, const char *varname)
{
  ppvar_t *v = pp_findvar(procpar, varname);
  if (!v) return UVAR_NOEXIST;
  return v->subtype;
}

enum ppvar_treemask pp_get_tree(const procpar_t *procpar, const char *varname)
{
  ppvar_t *v = pp_findvar(procpar, varname);
  if (!v) return TREE_NONE;
  return v->tree;
}

const char * pp_get_typename(enum ppvar_type type)
{
  switch (type) {
  case VAR_NOEXIST: return "NOEXIST";
  case VAR_INT: return "INT";
  case VAR_REAL: return "REAL";
  case VAR_STRING: return "STRING";
  default: return "UNKNOWN";
  }
}

const char * pp_get_subtypename(enum ppvar_subtype subtype)
{
  switch (subtype) {
  case UVAR_NOEXIST: return "NOEXIST";
  case UVAR_REAL: return "REAL";
  case UVAR_STRING: return "STRING";
  case UVAR_DELAY: return "DELAY";
  case UVAR_FLAG: return "FLAG";
  case UVAR_FREQ: return "FREQ";
  case UVAR_PULSE: return "PULSE";
  case UVAR_INT: return "INT";
  default: return "UNKNOWN";
  }
}

const char * pp_get_treename(enum ppvar_treemask tree)
{
  switch (tree) {
  case TREE_NONE: return "NONE";
  case TREE_GLOBAL: return "GLOBAL";
  case TREE_CURRENT : return "CURRENT";
  case TREE_PROC: return "PROC";
  case TREE_TMP: return "TMP";
  case TREE_SYSGLOBAL: return "SYSGLOBAL";
  case TREE_USER: return "USER";
  case TREE_ANY:
  default: return "UNKNOWN";
  }
}

/*! \brief	inquire array size for a variable
 */
int pp_get_count(procpar_t *procpar, const char *varname)
{
  ppvar_t *v = pp_findvar(procpar, varname);
  if (!v) return -1;
  return v->nvals;
}

/*! \brief	inquire array size for a variable's enums
 */
int pp_get_enumcount(procpar_t *procpar, const char *varname)
{
  ppvar_t *v = pp_findvar(procpar, varname);
  if (!v) return -1;
  return v->nenums;
}

/* \brief       return the parmater "varname"'s min field */
double pp_get_min(procpar_t *procpar, const char *varname)
{
  ppvar_t *v = pp_findvar(procpar, varname);
  if (!v) return -1;
  return v->min;
}

/* \brief       return the parmater "varname"'s max field */
double pp_get_max(procpar_t *procpar, const char *varname)
{
  ppvar_t *v = pp_findvar(procpar, varname);
  if (!v) return -1;
  return v->max;
}

/* \brief       return the parmater "varname"'s step field */
double pp_get_step(procpar_t *procpar, const char *varname)
{
  ppvar_t *v = pp_findvar(procpar, varname);
  if (!v) return -1;
  return v->step;
}

/* \brief       return the parmater "varname"'s group field */
int pp_get_group(procpar_t *procpar, const char *varname)
{
  ppvar_t *v = pp_findvar(procpar, varname);
  if (!v) return -1;
  return v->group;
}

/* \brief       return the parmater "varname"'s displaygroup field */
int pp_get_displaygroup(procpar_t *procpar, const char *varname)
{
  ppvar_t *v = pp_findvar(procpar, varname);
  if (!v) return -1;
  return v->display_group;
}


/*! \brief	get value of variable, if real. */
int pp_get_real(const procpar_t *procpar, const char *varname, int index, double *VAL)
{
  ppvar_t *v = pp_findvar(procpar, varname);
  if (!v) {
    //fprintf(stderr, "pp_get_real() accessing non-existant variable: '%s'[%d]\n", varname, index);
    *VAL = 0.0;
    return -1;
  }
  if (v->type != VAR_REAL && v->type != VAR_INT) {
    fprintf(stderr, "pp_get_real() ERROR: reading real from str '%s' (%d)\n", 
            varname, v->type);
    *VAL = 0.0;
    return -1;
  }
  assert(index >= 0 && index < v->nvals);
  *VAL = v->realvalues[index];
  return 0;
}

/*! \brief	get value of variable, if string.  return NULL if error */
int pp_get_string(const procpar_t *procpar, const char *varname, int index, const char ** STR)
{
  ppvar_t *v = pp_findvar(procpar, varname);
  if (!v) {
    //fprintf(stderr, "pp_get_string() accessing non-existant variable: '%s'[%d]\n", varname, index);
    *STR = NULL;
    return -1;
  }
  assert(v->type == VAR_STRING);
  assert(index >= 0 && index < v->nvals);
  *STR = v->strvalues[index];
  return 0;
}

/*! \brief	set value of variable, if real. */
int pp_set_real(procpar_t *procpar, const char *varname, int index, double value)
{
  ppvar_t *v = pp_findvar(procpar, varname);
  if (!v) {
    fprintf(stderr, "pp_set_real() accessing non-existant variable: '%s'[%d]\n", varname, index);
    return 0;
  }
  assert(v->type == VAR_REAL);
  assert(index >= 0);

  /* expand if needed */
  while (v->nvals <= index) {
    double * oldp = v->realvalues;
    v->nvals++;
    v->realvalues = realloc(v->realvalues, v->nvals * sizeof(double *));
    if (!v->realvalues) {
      fprintf(stderr, "pp_set_real() : realloc failed\n");
      v->realvalues = oldp;
      return -1;
    }
    v->realvalues[v->nvals - 1] = 0.;
  }

  v->realvalues[index] = value;
  return 0;
}

/*! \brief	set value of variable, if string.  return -1 if error */
int pp_set_string(procpar_t *procpar, const char *varname, int index, const char *value)
{
  ppvar_t *v = pp_findvar(procpar, varname);
  if (!v) {
    fprintf(stderr, "accessing non-existant variable: '%s'\n", varname);
    return -1;
  }
  assert(v->type == VAR_STRING);
  assert(index >= 0);

  /* expand if needed */
  while (v->nvals <= index) {
    char ** oldp = v->strvalues;
    v->nvals++;
    v->strvalues = realloc(v->strvalues, v->nvals * sizeof(char *));
    if (!v->strvalues) {
      fprintf(stderr, "pp_set_string() : realloc failed\n");
      v->strvalues = oldp;
      return -1;
    }
    v->strvalues[v->nvals - 1] = strdup(""); // very small potential memory leak here.
  }

  free(v->strvalues[index]);
  v->strvalues[index] = strdup(value);
  return 0;
}

/*! \brief      set an array of doubles into a "psuedo" variable
 *
 * \param[in]   procpar         procpar object to read from
 * \param[in]   varname         "psuedo" variable name
 * \param[out]  size            number of values to store
 * \param[out]  values          pointer to the values to store
 *
 * \return      0 for success, -1 for error
 */
int pp_set_realarray(procpar_t *procpar, const char *varname, int size, double *values)
{
  /* this stores a header variable called 'varname_header' with two
   * values in it -- 0=size, 1=number of splits, and then
   * splits the values into 'splits' different procpar variables
   * named varname_X, where X = {0,...,splits-1} */
  char varname2[MAX_STR];
  int splits;
  int splitsize;

  snprintf(varname2, sizeof(varname2), "%s_header", varname);
  splits = size /= (MAX_VARIABLE_CHARS / 18); // 18 ~ max chars in double string?
  splitsize = MAX_VARIABLE_CHARS / 18;
  splits = size / splitsize; // 18 ~ max chars in double string?

  if (pp_set_real(procpar, varname2, 0, size) ||
      pp_set_real(procpar, varname2, 1, splits)) {
    return -1;
  }
  int i, j;
  for (i = 0; i < splits; i++) {
    snprintf(varname2, sizeof(varname2), "%s_%d", varname, i);
    for (j = 0; j < splitsize; j++) {
      pp_set_real(procpar, varname2, j, values[i*splitsize + j]);
    }
  }
  return 0;
}

/*! \brief      set an array of doubles into a "psuedo" variable
 *
 * \param[in]   procpar         procpar object to read from
 * \param[in]   varname         "psuedo" variable name
 * \param[out]  size            receives the number of values
 * \param[out]  values          receives a malloc'd pointer to size doubles
 *
 * \return      0 for success, -1 for error
 */
int pp_get_realarray(procpar_t *procpar, const char *varname, int *size, double ** values)
{
  /* this stores a header variable called 'varname_header' with two
   * values in it -- 0=size, 1=number of splits, and then
   * splits the values into 'splits' different procpar variables
   * named varname_X, where X = {0,...,splits-1} */
  char varname2[MAX_STR];
  int splits;
  int splitsize;

  snprintf(varname2, sizeof(varname2), "%s_header", varname);

  if (VAR_REAL != pp_get_type(procpar, varname2))
    return -1;
  double tmp1;
  pp_get_real(procpar, varname2, 0, &tmp1);
  *size = tmp1;
  pp_get_real(procpar, varname2, 1, &tmp1);
  splits = tmp1;
  *values = malloc(sizeof(double) * *size);
  if (!*values)
    return -1;
  int curidx = 0;
  int i, j;
  for (i = 0; i < splits; i++) {
    snprintf(varname2, sizeof(varname2), "%s_%d", varname, i);
    if (VAR_REAL != pp_get_type(procpar, varname2))
      return -1;
    splitsize = pp_get_count(procpar, varname2);
    for (j = 0; j < splitsize && curidx < *size; j++) {
      pp_get_real(procpar, varname2, j, &(*values)[curidx]);
      curidx++;
    }
  }
  return 0;
}

/*! \brief      parse an assignment string and apply to the given procpar
 *
 * \return      0 for success, non-zero on error
 */
int pp_parse_assign(procpar_t *procpar, const char *string, enum ppvar_treemask tree)
{
  int err;
  regex_t preg;
  regmatch_t pmatch[10];
  const char *expr = "^([A-Za-z_][A-Za-z0-9_]*)(\\[([0-9]+)\\])?="
    "(([-+]?[0-9]+\\.?[0-9]*e?[-+]?[0-9]*)|(\"?([[:print:]]*)\"?))$";

  /*
   * run the string through the regex parser
   */

  //printf("regex: '%s' string: '%s'\n", expr, string);

  errno = 0;
  if ((err = regcomp(&preg, expr, REG_NEWLINE | REG_EXTENDED)) ||
      (err = regexec(&preg, string, sizeof(pmatch) / sizeof(pmatch[0]), pmatch, 0)))
  {
    char errbuf[2048];
    regerror(err, &preg, errbuf, sizeof(errbuf));
    fprintf(stderr, "pp_parse_assign: %s (%s ~= %s)\n", errbuf, expr, string);
    return -1;
  }

  /*
   * extract the meaningful bits of the string
   */
  char *varname;

  /* these match up to () positions in the expr */
  const int M_VARNAME = 1;
  const int M_INDEX = 3;
  const int M_REAL_VALUE = 5;
  const int M_STR_VALUE = 7;

#define MATCH_FOUND(idx) (pmatch[idx].rm_so != -1)
#define MATCH_DUP(idx) strndup(&string[pmatch[idx].rm_so], pmatch[idx].rm_eo - pmatch[idx].rm_so)

  /*
   * get the variable name
   */
  if (!MATCH_FOUND(M_VARNAME)) {
    fprintf(stderr, "no varname found in '%s'\n", string);
    return -1;
  }
  varname = MATCH_DUP(M_VARNAME);
  //fprintf(stderr, "str: '%s'\n", string);
  //fprintf(stderr, "varname rm_so:%ld  rm_eo:%ld  '%s'\n",
  //  (long int)pmatch[M_VARNAME].rm_so,
  //  (long int)pmatch[M_VARNAME].rm_eo, varname);

  /*
   * parse the index, default to 0 if not given
   */
  int index = 0;
  if (MATCH_FOUND(M_INDEX)) {
    char *indexstr;
    indexstr = MATCH_DUP(M_INDEX);
    printf("index: '%s'\n", indexstr);
    index = strtol_x(indexstr, 10, &err);
    if (err) {
      fprintf(stderr, "pp_parse_assign(): bad index: %s\n", indexstr);
      return -1;
    }
    free(indexstr);
  }

  enum ppvar_type valtype;
  enum ppvar_subtype valsubtype;
  char *valstr = NULL;
  double realval;
  int ret = 0;

  /*
   * parse the value to a {string, real, array}
   */
  if (MATCH_FOUND(M_STR_VALUE)) {
    valtype = VAR_STRING;
    valsubtype = UVAR_STRING;
    valstr = MATCH_DUP(M_STR_VALUE);
  }

  if (MATCH_FOUND(M_REAL_VALUE)) {
    valtype = VAR_REAL;
    valsubtype = UVAR_REAL;
    valstr = MATCH_DUP(M_REAL_VALUE);
    if (1 != sscanf(valstr, "%lg", &realval)) {
      fprintf(stderr, "pp_parse_assign(): real variable misformatted: %s\n", valstr);
      ret = -1;
      goto out;
    }
  }

  /*
   * do the assignment
   */
  enum ppvar_type type = pp_get_type(procpar, varname);
  //fprintf(stderr, "varname: '%s' value: '%s' type:%d\n", varname, valstr, valtype);

  /* doesn't exist, create a new variable in the procpar */
  if (type == VAR_NOEXIST) {
    ppvar_t *newvar = pp_newvar(varname, valtype, valsubtype, tree);
    pp_addvar(procpar, newvar);
    type = valtype;
  }

  switch (type) {
  case VAR_STRING:
    if (valtype != type) {
      fprintf(stderr, "wrong value type for variable: '%s' '%s'\n",
              varname, valstr);
      ret = -1;
      break;
    }
    ret = pp_set_string(procpar, varname, index, valstr);
    break;
  case VAR_INT:
  case VAR_REAL:
    if (valtype != type) {
      fprintf(stderr, "wrong value type for variable: '%s' '%s'\n",
              varname, valstr);
      ret = -1;
      break;
    }
    ret = pp_set_real(procpar, varname, index, realval);
    break;
  case VAR_NOEXIST: // cant happen
    break;
  }
out:
  free(varname);
  free(valstr);
  regfree(&preg);
  return ret;
}

#ifdef TEST

#include "cargs.h"

static argspec_t opts[] = {
  {'p', "procpar",      "procpar file to read", ARG_STRING},
  {'o', "outname",      "save modified procpar in this file", ARG_STRING},
  {'i', "varname",      "print a single variable from the procpar", ARG_STRING},
  {'v', "set_variable", "set variables, ie: -v np=512 -v sw=100000", ARG_STRING },
};

int main(int argc, char *argv[])
{
  char *ppfilename = NULL;
  char *outname = NULL;
  char *var = NULL;

  procpar_t *pp = pp_new();

  int i;
  argspec_t *opt = NULL;
  for (i = 0; i < argc; i++) {
    opt = cgetopt(argc, argv, NARRAY(opts), opts);
    if (!opt) {
      break;
    }
    switch (opt->flag) {
    case 'p': //
      ppfilename = opt->val.s;
      pp = pp_open(ppfilename, pp, TREE_CURRENT);
      if (!pp) {
        fprintf(stderr, "unable to open procpar: '%s'\n", ppfilename);
        return -1;
      }
      break;
    case 'o': //
      outname = opt->val.s;
      break;
    case 'i': //
      var = opt->val.s;
      break;
    case 'v':
      pp_parse_assign(pp, opt->val.s, TREE_CURRENT);
      break;
    }
  }

  if (!pp) {
    cerror(NARRAY(opts), opts, argv);
  }

  if (var)
    pp_print(pp, var);

  if (outname) {
    pp_write(outname, pp, TREE_ANY);
  }

  pp_close(pp);
  return 0;
}

#endif
