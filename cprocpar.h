/*
 * Copyright 2010-2014 Regents of the University of Minnesota.
 * This file subject to terms of Creative Commons CC BY-NC 4.0, see LICENSE.html
 *
 * author(s): Michael Tesch (tesch1@gmail.com)
 *
 */
#ifndef CPROCPAR_H
#define CPROCPAR_H

#define MAX_VARIABLE_CHARS (16384) /* estimate of vnmrj's character limitation */

#ifdef __cplusplus
extern "C" {
#endif

struct procpar_;
typedef struct procpar_ procpar_t;

/*! any particular variable's type can be:
 */
enum ppvar_type { VAR_NOEXIST=-1, VAR_INT=0, VAR_REAL=1, VAR_STRING=2 };
enum ppvar_subtype { UVAR_NOEXIST=-1, UVAR_REAL=1, UVAR_STRING=2, 
                     UVAR_DELAY=3, UVAR_FLAG=4, UVAR_FREQ=5, UVAR_PULSE=6,
                     UVAR_INT=7};
enum ppvar_treemask { TREE_NONE=        0x80000000,
                      TREE_GLOBAL=      0x00000001,
                      TREE_CURRENT=     0x00000002,
                      TREE_PROC=        0x00000004,
                      TREE_TMP=         0x00000008,
                      TREE_SYSGLOBAL=   0x00000010,
                      TREE_USER=        0x00000020,
                      TREE_ANY=         0x0000003f};

procpar_t * pp_new();
procpar_t * pp_open(const char *filename, procpar_t * tomerge, enum ppvar_treemask tree);
procpar_t * pp_merge(procpar_t * dest, procpar_t * old);
void pp_close(procpar_t *procpar);
int pp_write(const char *filename, procpar_t *procpar, enum ppvar_treemask tree);
int pp_delete(procpar_t *procpar, const char * varname);
int pp_size(procpar_t *procpar);
const char * pp_get_next(procpar_t *procpar);

enum ppvar_type pp_get_type(const procpar_t *procpar, const char *varname);
enum ppvar_subtype pp_get_subtype(const procpar_t *procpar, const char *varname);
enum ppvar_treemask pp_get_tree(const procpar_t *procpar, const char *varname);
const char * pp_get_typename(enum ppvar_type type);
const char * pp_get_subtypename(enum ppvar_subtype type);
const char * pp_get_treename(enum ppvar_treemask type);
int pp_get_count(procpar_t *procpar, const char *varname);
int pp_get_enumcount(procpar_t *procpar, const char *varname);
int pp_get_count(procpar_t *procpar, const char *varname);
double pp_get_min(procpar_t *procpar, const char *varname);
double pp_get_max(procpar_t *procpar, const char *varname);
double pp_get_step(procpar_t *procpar, const char *varname);
int pp_get_group(procpar_t *procpar, const char *varname);
int pp_get_displaygroup(procpar_t *procpar, const char *varname);

int pp_get_real(const procpar_t *procpar, const char *varname, int index, double *VAL);
int pp_get_string(const procpar_t *procpar, const char *varname, int index, const char ** STR);
int pp_set_real(procpar_t *procpar, const char *varname, int index, double value);
int pp_set_string(procpar_t *procpar, const char *varname, int index, const char *value);
int pp_set_realarray(procpar_t *procpar, const char *varname, int size, double *values);
int pp_get_realarray(procpar_t *procpar, const char *varname, int *size, double ** values);

int pp_parse_assign(procpar_t *procpar, const char *string, enum ppvar_treemask tree);

#ifdef __cplusplus
}
#endif
#endif /* CPROCPAR_H */
