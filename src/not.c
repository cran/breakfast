#include "not.h"

SEXP contrasts_t_to_dataframe(contrasts_t *contrasts){
  
  SEXP lst, start, end, arg_max, max, lstnames, len;
  
  PROTECT(lst = allocVector(VECSXP, 5));
  PROTECT(start = allocVector(INTSXP, (*contrasts).n_intervals));
  PROTECT(end = allocVector(INTSXP, (*contrasts).n_intervals));
  PROTECT(arg_max = allocVector(INTSXP, (*contrasts).n_intervals));
  PROTECT(max = allocVector(REALSXP, (*contrasts).n_intervals));
  PROTECT(len = allocVector(INTSXP, (*contrasts).n_intervals));
  PROTECT(lstnames = allocVector(STRSXP, 5));
  
  
  SET_VECTOR_ELT(lst, 0, start);
  SET_VECTOR_ELT(lst, 1, end);
  SET_VECTOR_ELT(lst, 2, len);
  SET_VECTOR_ELT(lst, 3, arg_max);
  SET_VECTOR_ELT(lst, 4, max);
  
  SET_STRING_ELT(lstnames, 0, mkChar("start")); 
  SET_STRING_ELT(lstnames, 1, mkChar("end")); 
  SET_STRING_ELT(lstnames, 2, mkChar("length")); 
  SET_STRING_ELT(lstnames, 3, mkChar("arg max")); 
  SET_STRING_ELT(lstnames, 4, mkChar("max contrast")); 
  
  setAttrib(lst,R_NamesSymbol,lstnames); 
  
  memcpy(INTEGER(start), (*contrasts).start, (*contrasts).n_intervals * sizeof(int));
  memcpy(INTEGER(end), (*contrasts).end, (*contrasts).n_intervals * sizeof(int));
  memcpy(INTEGER(len), (*contrasts).length, (*contrasts).n_intervals * sizeof(int));
  memcpy(INTEGER(arg_max), (*contrasts).arg_max, (*contrasts).n_intervals * sizeof(int));
  memcpy(REAL(max), (*contrasts).max, (*contrasts).n_intervals * sizeof(double));
  
  SEXP df;
  df = PROTECT(lang2(install("data.frame"), lst)); 
  SEXP res = PROTECT(eval(df, R_GlobalEnv)); 
  
  UNPROTECT(9);
  
  return res;
}


SEXP solution_path_t_to_list(solution_path_t *solution_path){
  
  
  SEXP res, resnames, tmp_cpt, tmp_index, th, n_cpt, locations, index;
  int i;
  int n_protected = 0;
  int n_th = (*solution_path).n_th;
  
  PROTECT(locations = allocVector(VECSXP, n_th));
  PROTECT(index = allocVector(VECSXP, n_th));  
  PROTECT(th = allocVector(REALSXP, n_th));
  PROTECT(n_cpt = allocVector(INTSXP, n_th));
  
  n_protected += 4;
  
  //copy the estimated change-points
  
  double *p_th = REAL(th);
  int *p_n_cpt = INTEGER(n_cpt);
  
  for(i=0; i<n_th; i++){
    
    p_th[n_th-i-1] = (*solution_path).cpts[i].min_max - DBL_EPSILON;
    p_n_cpt[n_th-i-1] = (*solution_path).cpts[i].n_cpt;
    
    PROTECT(tmp_cpt = allocVector(INTSXP, (*solution_path).cpts[i].n_cpt));
    memcpy(INTEGER(tmp_cpt), (*solution_path).cpts[i].cpt, (*solution_path).cpts[i].n_cpt * sizeof(int));
    SET_VECTOR_ELT(locations, n_th-i-1, tmp_cpt);
    
    PROTECT(tmp_index = allocVector(INTSXP, (*solution_path).cpts[i].n_cpt));
    memcpy(INTEGER(tmp_index), (*solution_path).cpts[i].index, (*solution_path).cpts[i].n_cpt * sizeof(int));
    SET_VECTOR_ELT(index, n_th-i-1, tmp_index);
    
    UNPROTECT(2);
  }
  
  //create the list with results
  PROTECT(res = allocVector(VECSXP, 4));
  n_protected++;
  
  SET_VECTOR_ELT(res, 0, th);
  SET_VECTOR_ELT(res, 1, locations);
  SET_VECTOR_ELT(res, 2, index);  
  SET_VECTOR_ELT(res, 3, n_cpt);
  PROTECT(resnames = allocVector(STRSXP, 4));
  n_protected++;
  
  SET_STRING_ELT(resnames,0,mkChar("th")); 
  SET_STRING_ELT(resnames,1,mkChar("cpt")); 
  SET_STRING_ELT(resnames,2,mkChar("index"));
  SET_STRING_ELT(resnames,3,mkChar("n.cpt"));
  
  setAttrib(res, R_NamesSymbol, resnames); 
  
  //Clean up
  
  UNPROTECT(n_protected);
  
  //return the outcome
  return res;
}



SEXP not_r_wrapper(SEXP x, SEXP intervals, SEXP method, SEXP contrast_type, SEXP parallel, SEXP augmented){
  
  SEXP intervals_dim;
  int v_n_protected = 0;
  
  PROTECT(intervals_dim = getAttrib(intervals, R_DimSymbol));
  v_n_protected++;
  
  int v_n_obs = length(x);
  int v_n_intervals = INTEGER(intervals_dim)[0];

  contrasts_t *p_contrasts;
  double *p_x = REAL(x);
  int *p_intervals = (int *) INTEGER(intervals);
  
  int v_parallel = INTEGER(parallel)[0];
  int v_method = INTEGER(method)[0];
  int v_contrast_type = INTEGER(contrast_type)[0];
  int v_augmented = INTEGER(augmented)[0];
  int v_min_dist;
  
  
  // Find CUSUMS over randomly drawn intervals
  eval_contrast_fun_t eval_contrast_fun;
  
  switch(v_contrast_type){
    case 1: 
      eval_contrast_fun = &slope_contrast;
      v_min_dist = 2;
      break;
    case 3: 
      eval_contrast_fun = &intercept_slope_and_quadratic_contrast;
      v_min_dist = 3;
      break;
    case 2:
      eval_contrast_fun = &intercept_and_slope_contrast;
      v_min_dist = 2;
      break;
    case 4: 
      eval_contrast_fun = &intercept_and_volatility_contrast;
      v_min_dist = 2;
      break;
    case 5: 
      eval_contrast_fun = &intercept_signs_contrast;
      v_min_dist = 1;
    break;
    default:
      eval_contrast_fun = &intercept_contrast;
      v_min_dist = 1;
  }
  
  p_contrasts = eval_contrasts(p_x, v_n_obs, p_intervals, v_n_intervals, eval_contrast_fun, v_parallel);
  
  
  // Convert p_contrasts to an R object
  SEXP contrasts;
  
  PROTECT(contrasts = contrasts_t_to_dataframe(p_contrasts));
  v_n_protected++;
  
  
  // Depending on the algorithm, either sort the contrasts according to either their values or interval lengths
  double * p_tmp = Calloc((*p_contrasts).n_intervals, double);
  
  // v_method = 0 -> not
  // v_method > 1 -> max i.e. wbs
  
  if(v_method == 0){
    for(int i=0; i<(*p_contrasts).n_intervals; i++) p_tmp[i] =  (double) ((*p_contrasts).length[i]);
    // Sort the contrasts from the shortest intervals to the largest one
    rsort_with_index(p_tmp, (int *)((*p_contrasts).index), (*p_contrasts).n_intervals);
  }else{
    // Sort the contrasts from the largest val to the smallest one
    for(int i=0; i<(*p_contrasts).n_intervals; i++) p_tmp[i] =  (double) ((*p_contrasts).max[i]);
    revsort(p_tmp, (int *)((*p_contrasts).index), (*p_contrasts).n_intervals);
  }
  
  //for(int i=0; i<(*p_contrasts).n_intervals; i++) Rprintf("%d ", (*p_contrasts).index[i]+1);
  //Rprintf("\n");
  
  Free(p_tmp);
  
  solution_path_t *p_solution_path; 
  
  // Find change-points for all possible thresholds
  if(v_augmented == 0) p_solution_path = solution_path(p_contrasts, NULL, v_min_dist);
  //// suppress the augmented option here for now
  // else p_solution_path = solution_path(p_contrasts, eval_contrast_fun, v_min_dist);
  else p_solution_path = solution_path(p_contrasts, NULL, v_min_dist);
  
  
  // Convert  p_solution_path to an R object
  SEXP sol_path;
  PROTECT(sol_path = solution_path_t_to_list(p_solution_path));
  v_n_protected++;
  
  // Create list with solution path and results
  SEXP res, resnames;
  
  PROTECT(res = allocVector(VECSXP,2));
  PROTECT(resnames = allocVector(STRSXP,2));
  v_n_protected += 2;
  
  SET_VECTOR_ELT(res, 0, contrasts);
  SET_VECTOR_ELT(res, 1, sol_path);
  
  SET_STRING_ELT(resnames, 0, mkChar("contrasts"));
  SET_STRING_ELT(resnames, 1, mkChar("solution.path"));
  
  
  setAttrib(res,R_NamesSymbol,resnames); 
  
  
  //Clean up...
  destroy_solution_path(&p_solution_path);
  destroy_contrasts(&p_contrasts);
  UNPROTECT(v_n_protected);

  //return the result
  
  return res;
  
  
}


