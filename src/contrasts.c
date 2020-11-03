#include "contrasts.h"

//functions for change-point detection in piecewise constant signal + noise
void alloc_contrasts(contrasts_t **contrasts, int n_intervals, double *x, int n_obs){
  
  (*contrasts) = Calloc(1, contrasts_t);
  (**contrasts).index = Calloc(n_intervals, int);
  (**contrasts).start = Calloc(n_intervals, int);
  (**contrasts).end = Calloc(n_intervals, int);
  (**contrasts).length = Calloc(n_intervals, int);
  (**contrasts).arg_max = Calloc(n_intervals, int);
  (**contrasts).max = Calloc(n_intervals, double);
  (**contrasts).x = Calloc(n_obs, double);
  memcpy((**contrasts).x, x, n_obs * sizeof(double));
  (**contrasts).n_obs = n_obs;
  (**contrasts).n_intervals = n_intervals;
  
}

void destroy_contrasts(contrasts_t **contrasts){
  
  if( (*contrasts) != NULL){
    
    Free((**contrasts).index);
    Free((**contrasts).start);
    Free((**contrasts).end);
    Free((**contrasts).index);
    Free((**contrasts).arg_max);
    Free((**contrasts).max);
    Free((**contrasts).x);
    Free(*contrasts);
    
    (*contrasts) = NULL;
    
  }
  
}


max_contrast_t intercept_contrast(double *x, int n_obs){
  
  max_contrast_t max_contrast;
  max_contrast.arg_max = 0;
  max_contrast.max = 0;
  
  if(n_obs > 1){

    int i;
    
    double sum_left=0, sum_right=0, const_coef;
    
    sum_left = x[0];
    for(i=1; i < n_obs; i++) sum_right += x[i];
    
    double n_dbl = (double) n_obs;
    double c=1, n_min_c=n_dbl-1;
    
    const_coef = n_min_c * sum_left - c * sum_right;
    const_coef *= const_coef;
    const_coef *= (1.0/(n_dbl * n_min_c * c));
    
    max_contrast.max = const_coef;

    
    for(i=1; i<(n_obs-1); i++){
      
      sum_left += x[i];
      sum_right -= x[i];
      
      c = i+1;
      n_min_c = n_dbl - c;
      
      const_coef = n_min_c * sum_left - c * sum_right;
      const_coef *= const_coef;
      const_coef *= (1.0 / (n_dbl * n_min_c * c));
      

      
      if(const_coef > max_contrast.max) {
        max_contrast.max = const_coef;
        max_contrast.arg_max = i;
      }
      
    
    }
    
    max_contrast.max = sqrt(max_contrast.max);
    

  }
  
  return max_contrast;
  
}

max_contrast_t intercept_signs_contrast(double *x, int n_obs){
  
  max_contrast_t max_contrast;
  max_contrast.arg_max = 0;
  max_contrast.max = 0;
  
  if(n_obs > 3){
    
    int i;
    
    //find the mean in this interval
    double mean = 0.0;
    
    for(i=0; i < n_obs; i++) mean += x[i];
    mean /= n_obs;
    
    double sum_left=0, sum_right=0, const_coef;
    
    sum_left = sign(x[0]-mean);
    for(i=1; i < n_obs; i++) sum_right += sign(x[i]-mean);
    
    double n_dbl = (double) n_obs;
    double c=1, n_min_c=n_dbl-1;
    
    const_coef = n_min_c * sum_left - c * sum_right;
    const_coef *= const_coef;
    const_coef *= (1.0/(n_dbl * n_min_c * c));
    
    max_contrast.max = const_coef;
    
    
    for(i=1; i<(n_obs-1); i++){
      
      sum_left += sign(x[i]-mean);
      sum_right -= sign(x[i]-mean);
      
      c = i+1;
      n_min_c = n_dbl - c;
      
      const_coef = n_min_c * sum_left - c * sum_right;
      const_coef *= const_coef;
      const_coef *= (1.0 / (n_dbl * n_min_c * c));
      
      
      
      if(const_coef >= max_contrast.max) {
        max_contrast.max = const_coef;
        max_contrast.arg_max = i;
      }
      
      
    }
    
    max_contrast.max = sqrt(max_contrast.max);
    
    
  }
  
  return max_contrast;
  
}

max_contrast_t slope_contrast(double *x, int n_obs){
  
  max_contrast_t max_contrast;
  max_contrast.arg_max = 0;
  max_contrast.max = 0.0;
  
  if(n_obs > 3){
    
    
    int i,j;
    
    double lin_coef;
    double n_dbl = (double) n_obs;
    double c;
    
    double *sum_left = Calloc(n_obs, double);
    double *sum_right = Calloc(n_obs, double);
    double *lin_sum_left = Calloc(n_obs, double);
    double *lin_sum_right = Calloc(n_obs, double);
    
    sum_left[0] = x[0];
    lin_sum_left[0] = x[0];
    
    j = n_obs-1;
    sum_right[j] = x[j];
    lin_sum_right[j] = x[j] * n_obs;
    
    
    for(i=1; i<n_obs; i++){
      
      sum_left[i] = sum_left[i-1] + x[i];
      lin_sum_left[i] = lin_sum_left[i-1] + (i+1) * x[i];
      
      j = n_obs-i-1;
      sum_right[j] = sum_right[j+1] + x[j];
      lin_sum_right[j] = lin_sum_right[j+1] + (j+1) * x[j];
      
    }

    max_contrast.max = 0;
    max_contrast.arg_max = 0;
    
    double constant_factor = 6.0/((n_dbl - 1)*(n_dbl)*(n_dbl+1));
    double n_min_c, tmp1, tmp2;
      
    for(i=2; i<(n_obs-2); i++){
      
      c = (double)(i+1);
      n_min_c = n_dbl - c;
      
      tmp1 = (n_min_c+1) * n_min_c;
      tmp2 = c * (c-1);
      
      lin_coef = tmp1 * (lin_sum_left[i]  * (n_dbl + 2*c - 1) - sum_left[i] * (c*n_dbl + c + 1));
      lin_coef -= tmp2 * (lin_sum_right[i+1] *(3*n_dbl - 2*c + 1)- sum_right[i+1]*(2*n_dbl - c*n_dbl + 2*n_dbl*n_dbl - c));
      lin_coef *= lin_coef;
      lin_coef *= constant_factor * (1.0 / (tmp1 * tmp2)) * (1.0 / ((1-2*c*c+2*c*n_dbl+2*c-n_dbl)));
   
      if(lin_coef >= max_contrast.max){
        max_contrast.max = lin_coef;
        max_contrast.arg_max = i;
      }
      
    }
    
    max_contrast.max = sqrt(max_contrast.max);
    
    Free(sum_left);
    Free(sum_right);
    Free(lin_sum_left);
    Free(lin_sum_right);
    
  }  
  
  return max_contrast;
  
}


max_contrast_t intercept_and_slope_contrast(double *x, int n_obs){
  
  max_contrast_t max_contrast;
  max_contrast.arg_max = 0;
  max_contrast.max = 0;

  if(n_obs > 5){
    
    int i,j;
    
    double contrast, const_coef, lin_coef, lin_coef_left, lin_coef_right;
    double n_dbl = (double) n_obs;
    double c, n_min_c;
    
    max_contrast.max = 0;
    max_contrast.arg_max = 0;
    
    double *sum_left = Calloc(n_obs, double);
    double *sum_right = Calloc(n_obs, double);
    double *lin_sum_left = Calloc(n_obs, double);
    double *lin_sum_right = Calloc(n_obs, double);
    
    sum_left[0] = x[0];
    lin_sum_left[0] = x[0];
    
    j = n_obs-1;
    sum_right[j] = x[j];
    lin_sum_right[j] = x[j] * n_obs;
    
    
    for(i=1; i<n_obs; i++){
      
      sum_left[i] = sum_left[i-1] + x[i];
      lin_sum_left[i] = lin_sum_left[i-1] + (i+1) * x[i];
      
      j = n_obs-i-1;
      sum_right[j] = sum_right[j+1] + x[j];
      lin_sum_right[j] = lin_sum_right[j+1] + (j+1) * x[j];
      
    }
    


    lin_coef =  (2.0 * lin_sum_left[n_obs-1] - (((n_dbl+1)) * sum_left[n_obs-1]));
    lin_coef *= lin_coef;
    lin_coef /= n_dbl *(n_dbl * n_dbl -1) / 3.0;


    for(i=1; i< (n_obs-2); i++){
      
      c = (double) (i+1);
      n_min_c = n_dbl - c;

      const_coef = n_min_c * sum_left[i] - c * sum_right[i+1];
      const_coef *= const_coef;
      const_coef /= n_dbl * n_min_c * c;
      
      
      lin_coef_left = (2.0 * lin_sum_left[i] - ((c+1) * sum_left[i]));
      lin_coef_left *= lin_coef_left;
      lin_coef_left /= c*(c*c-1) / 3;
      
      
      lin_coef_right =  (2.0 * lin_sum_right[i+1] - ((c+1+n_dbl) * sum_right[i+1]));
      lin_coef_right *= lin_coef_right;
      lin_coef_right /= n_min_c * (n_min_c*n_min_c-1) / 3.0;
      
      
      contrast = const_coef + lin_coef_left + lin_coef_right - lin_coef;
      
      if(contrast >= max_contrast.max){
        max_contrast.max = contrast;
        max_contrast.arg_max = i;
      }
      
    }
    
    max_contrast.max = sqrt(max_contrast.max);
      
    Free(sum_left);
    Free(sum_right);
    Free(lin_sum_left);
    Free(lin_sum_right);
    
  }  

  
  return max_contrast;
  
}

max_contrast_t intercept_slope_and_quadratic_contrast(double *x, int n_obs){
  
  max_contrast_t max_contrast;
  max_contrast.max = 0;
  max_contrast.arg_max = 0;
  
  if(n_obs > 7){
    
  
    int i,j;
    
    double contrast, const_coef, lin_coef, lin_coef_left, lin_coef_right, quad_coef, quad_coef_left, quad_coef_right;
    double n_dbl = (double) n_obs;
    double c, n_min_c;

    double *sum_left = Calloc(n_obs, double);
    double *sum_right = Calloc(n_obs, double);
    double *lin_sum_left = Calloc(n_obs, double);
    double *lin_sum_right = Calloc(n_obs, double);
    double *quad_sum_left = Calloc(n_obs, double);
    double *quad_sum_right = Calloc(n_obs, double);

    
    sum_left[0] = x[0];
    lin_sum_left[0] = x[0];
    quad_sum_left[0] = x[0];
    
    j = n_obs-1;
    sum_right[j] = x[j];
    lin_sum_right[j] = x[j] * n_obs;
    quad_sum_right[j] = x[j] * n_obs * n_obs;
    
    for(i=1; i<n_obs; i++){
      
      sum_left[i] = sum_left[i-1] + x[i];
      lin_sum_left[i] = lin_sum_left[i-1] + (i+1) * x[i];
      quad_sum_left[i] = quad_sum_left[i-1] + (i+1) * (i+1) * x[i];
      
      j = n_obs-i-1;
      sum_right[j] = sum_right[j+1] + x[j];
      lin_sum_right[j] = lin_sum_right[j+1] + (j+1) * x[j];
      quad_sum_right[j] = quad_sum_right[j+1] + (j+1) * (j+1) * x[j];
      
    }
    
    
    lin_coef =  (2.0 * lin_sum_left[n_obs-1] - (((n_dbl+1)) * sum_left[n_obs-1]));
    lin_coef *= lin_coef;
    lin_coef *= (3.0 / (n_dbl * (n_dbl*n_dbl-1)));
    
    
    
    quad_coef = (6 * quad_sum_left[n_obs-1] - (6*n_dbl +6 ) * lin_sum_left[n_obs-1]  +  (2 + 3*n_dbl + n_dbl*n_dbl) * sum_left[n_obs-1]);
    quad_coef *= quad_coef;
    quad_coef *= (5.0 / (n_dbl * (n_dbl*n_dbl-1) * (n_dbl*n_dbl-2)));;
    

    
    for(i=2; i<(n_obs-3); i++){
      
      c = (double) (i+1);
      n_min_c = n_dbl - c;
      
      const_coef = n_min_c * sum_left[i] - c * sum_right[i+1];
      const_coef *= const_coef;
      const_coef *= (1.0/(n_dbl * n_min_c * c));
      

      
      lin_coef_left =  (2.0 * lin_sum_left[i] - (((c+1)) * sum_left[i]));
      lin_coef_left *= lin_coef_left;
      lin_coef_left *= (3.0 / (c * (c*c-1)));
      
      lin_coef_right =  (2.0 * lin_sum_right[i+1] - (((c+1+n_dbl)) * sum_right[i+1] ));
      lin_coef_right *= lin_coef_right;
      lin_coef_right *= (3.0 / (n_min_c * (n_min_c*n_min_c-1)));
      
      quad_coef_left = (6 * quad_sum_left[i] - (6*c +6 ) * lin_sum_left[i]  +  (2 + 3*c + c*c) * sum_left[i]);
      quad_coef_left *= quad_coef;
      quad_coef_left *= (5.0 / (c * (c * c-1) * (c*c-2)));

      
      quad_coef_right = 6 * quad_sum_right[i+1] - 6 * (1+c+n_dbl) * lin_sum_right[i+1];
      quad_coef_right += (2 + c * (3 + c + 4 * n_dbl) + n_dbl * (n_dbl +3)) * sum_right[i+1];
      quad_coef_right *= quad_coef;
      quad_coef_right *=  (5.0 / (n_min_c * (n_min_c * n_min_c-1) * (n_min_c*n_min_c-2)));    
           
            
      contrast = const_coef + lin_coef_left + lin_coef_right - lin_coef + quad_coef_left + quad_coef_right - quad_coef;
      
      if(contrast >= max_contrast.max){
        max_contrast.max = contrast;
        max_contrast.arg_max = i;
      }

      
    }
    
    max_contrast.max = sqrt(max_contrast.max);
  
    Free(sum_left);
    Free(sum_right);
    Free(lin_sum_left);
    Free(lin_sum_right);
    Free(quad_sum_left);
    Free(quad_sum_right);

  }
  
  return max_contrast;
  
}

max_contrast_t intercept_and_volatility_contrast(double *x, int n_obs){
  
  max_contrast_t max_contrast;
  max_contrast.arg_max = 0;
  max_contrast.max = 0;
  
  if(n_obs >=6){
    
    int i,j;
    
    double contrast, var, var_left, n_log_var, var_right, tmp;
    double n_dbl = (double) n_obs;
    double c, n_min_c;
    double eps = sqrt(DBL_EPSILON);
    
    max_contrast.max = 0;
    max_contrast.arg_max = 0;
    
    double *sum_left = Calloc(n_obs, double);
    double *sum_right = Calloc(n_obs, double);
    double *sum_sq_left = Calloc(n_obs, double);
    double *sum_sq_right = Calloc(n_obs, double);
    
    sum_left[0] = x[0];
    sum_sq_left[0] = x[0] * x[0];
    
    j = n_obs-1;
    sum_right[j] = x[j];
    sum_sq_right[j] = x[j] * x[j];
    
    
    for(i=1; i<n_obs; i++){
      
      sum_left[i] = sum_left[i-1] + x[i];
      sum_sq_left[i] = sum_sq_left[i-1] +  x[i] * x[i];
      
      j = n_obs-i-1;
      sum_right[j] = sum_right[j+1] + x[j];
      sum_sq_right[j] = sum_sq_right[j+1] + x[j] * x[j];
      
    }
    
    
    tmp = sum_right[0]  / n_dbl;
    var = sum_sq_right[0] / n_dbl  - tmp * tmp;
    
    
    if (fabs(var) < eps){
      
      max_contrast.arg_max = n_obs / 2;
      
    } else{
      
      n_log_var = n_dbl * log(var);
        
      for(i=3; i < (n_obs-4); i++){
        
        c = (double)(i+1);
        n_min_c = n_dbl - c;
        
        tmp = sum_left[i]  / c;
        var_left = sum_sq_left[i] / c  - tmp * tmp;
        
        tmp = sum_right[i+1]  / n_min_c;
        var_right = sum_sq_right[i+1] / n_min_c  - tmp * tmp;
        
        if(fabs(var_left) < eps || fabs(var_right) < eps ) contrast = 0;
        else contrast = -2 * (c * log(var_left) + n_min_c * log(var_right) - n_log_var);
        
        
        if(contrast >= max_contrast.max){
          max_contrast.max = contrast;
          max_contrast.arg_max = i;
        }
        
      }
      
    }
    
    Free(sum_left);
    Free(sum_right);
    Free(sum_sq_left);
    Free(sum_sq_right);
    
  }  
  
  return max_contrast;
  
}

contrasts_t *eval_contrasts(double *x, int n_obs, int *intervals, int n_intervals,
                            eval_contrast_fun_t eval_contrast_fun, int parallel){
  contrasts_t *contrasts;
  alloc_contrasts(&contrasts, n_intervals, x, n_obs);
  int *start = intervals; 
  int *end = &intervals[n_intervals];
  
  
  //computation using one core only
  if(parallel == 0){
    
    int i;
    max_contrast_t max_contrast;
    int interval_length;
    
    
    for(i=0; i<n_intervals; i++){
      
      interval_length =  end[i]-start[i]+1;
      
      max_contrast = eval_contrast_fun(&x[start[i]-1], interval_length);
      
      (*contrasts).start[i] = start[i];
      (*contrasts).end[i] = end[i];
      (*contrasts).length[i] = interval_length;
      (*contrasts).max[i] = max_contrast.max;
      (*contrasts).arg_max[i] = max_contrast.arg_max +start[i];
      (*contrasts).index[i] = i;
      
    }
    
    
  }else{ //computations using all available cores
    
    #pragma omp parallel
    {
      int i;
      max_contrast_t max_contrast;
      int interval_length;
      
      #pragma omp for
      for(i=0; i<n_intervals; i++){
        
        interval_length =  end[i]-start[i]+1;
        
        max_contrast = eval_contrast_fun(&x[start[i]-1], interval_length);
        
        (*contrasts).start[i] = start[i];
        (*contrasts).end[i] = end[i];
        (*contrasts).length[i] = interval_length;
        (*contrasts).max[i] = max_contrast.max;
        (*contrasts).arg_max[i] = max_contrast.arg_max +start[i];
        (*contrasts).index[i] = i;
        
      }
        
      
    }
  }
  
  return contrasts;
}
