#ifndef NOT_H
#define NOT_H

#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>


#include "contrasts.h"
#include "changepoints_tree.h"

#define IDX(i,j,ld) ((((j)-1) * (ld))+((i)-1))

typedef struct ip_max{
  int arg_max;
  double max;
  double abs_max;
} ip_max_t;

typedef struct ips {
  int *index;
  int *s;
  int *e;
  int *cpt;
  double *max;
  double *abs_max;
  int M;
  int n;
} ips_t;


typedef struct notres{
  int *s;
  int *e;
  int *cpt;
  double *max;
  double *minth;
  int *scale;
  int M;
  int n;
} not_res_t;

SEXP not_r_wrapper(SEXP x, SEXP intervals, SEXP method, SEXP contrast_type, SEXP parallel, SEXP augmented);
SEXP solution_path_t_to_list(solution_path_t *solution_path);
SEXP contrasts_t_to_dataframe(contrasts_t *contrasts);


#endif
