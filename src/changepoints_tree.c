#include "changepoints_tree.h"

void build_tree(cpt_tree_node_t **node, cpt_tree_node_t **parent_node,  int start, int end, double th, contrasts_t *contrasts, eval_contrast_fun_t eval_contrast_fun){
  
  
 if(end-start>=1) {
  // first case - the node is empty
  if((*node) == NULL){
    
    // find the indices of the intervals such that they are greater than threshold and they are inside (s, e)
    // if there is no parent node - create index, if the parent node exists, use it's index
    
    int *index;
    int n_intervals;
    
    if((*parent_node) == NULL){
      
      n_intervals = (*contrasts).n_intervals;
      index =(*contrasts).index;
      
    }else{
      
      n_intervals = (**parent_node).n_intervals;
      index = (**parent_node).index;
      
    }
    
     
    int *new_index = Calloc(n_intervals, int);
    int i,j,new_n_intervals=0;
    
    for(i=0; i<n_intervals; i++){
      
      j = index[i];
      if( ((*contrasts).max[j] > th) &&  ((*contrasts).start[j] >= start) && ((*contrasts).end[j] <= end)) new_index[new_n_intervals++] = j;
      
    }
    
    if(eval_contrast_fun == NULL){
      
     //if nothing found - clean up, else-create new nodes
      if(new_n_intervals == 0) Free(new_index);
      else{
        
        new_index = Realloc(new_index, new_n_intervals, int);
        (*node) = Calloc(1, cpt_tree_node_t);
        
  
        //struct cpt_tree_node *left_node, *right_node;
        
        (**node).index = new_index;
        (**node).n_intervals = new_n_intervals;
        (**node).left_node = NULL;
        (**node).right_node = NULL;
        
        j = new_index[0];
        
        (**node).cpt = (*contrasts).arg_max[j];
        (**node).max = (*contrasts).max[j];
        
        //build trees for the children nodes
        build_tree(&((**node).left_node), node, start, (**node).cpt, th, contrasts, eval_contrast_fun);
        build_tree(&((**node).right_node), node, (**node).cpt+1, end, th, contrasts, eval_contrast_fun);
        
      }
      
    }
    
    // NOTE: to make it coherent with other parts of the package
    // the augmented option is suppressed here for now, which combines BS with WBS 
    // - i.e. the whole interval would always be considered there
    /* ------

    else  {
      // Modified wbs in case number of drawn intervals is small
      // i.e. Combining bs and wbs into one
      // so always consider the whole subinterval as well
      max_contrast_t contrast = eval_contrast_fun(&((*contrasts).x[start-1]), end-start+1);
      
      if(new_n_intervals == 0) Free(new_index);
      
      if((contrast.max > th) & (new_n_intervals==0)){

        (*node) = Calloc(1, cpt_tree_node_t);
        (**node).index = NULL;
        (**node).n_intervals = 0;
        (**node).left_node = NULL;
        (**node).right_node = NULL;
        
        (**node).cpt = contrast.arg_max+start;
        (**node).max = contrast.max;
        
        build_tree(&((**node).left_node), node, start, (**node).cpt, th, contrasts, eval_contrast_fun);
        build_tree(&((**node).right_node), node, (**node).cpt+1, end, th, contrasts, eval_contrast_fun);

        
      } else if(new_n_intervals > 0){
        
        
        new_index = Realloc(new_index, new_n_intervals, int);
        
        (*node) = Calloc(1, cpt_tree_node_t);
        (**node).index = new_index;
        (**node).n_intervals = new_n_intervals;
        (**node).left_node = NULL;
        (**node).right_node = NULL;
        
        j = new_index[0];

        if((*contrasts).max[j] > contrast.max) {
          (**node).cpt = (*contrasts).arg_max[j];  
          (**node).max = (*contrasts).max[j];
        }else{
          (**node).cpt = contrast.arg_max+start;;
          (**node).max = contrast.max;
        }
        
        build_tree(&((**node).left_node), node, start, (**node).cpt, th, contrasts, eval_contrast_fun);
        build_tree(&((**node).right_node), node, (**node).cpt+1, end, th, contrasts, eval_contrast_fun);
        
      }
    }
    ------*/
  } 

    
  else {
    //this when we just modify the tree applying a new threshold 
    
    //we check if we need to modify the tree
    if((**node).max <= th) {
      
      //the entire branch of the tree needs to be reconstructed, hence we destroy it
      destroy_tree(node);

      //we build this branch from scratch
      build_tree(node, parent_node, start, end, th, contrasts, eval_contrast_fun);
      
    } else {
      
      //there is no need to rebuild the tree - we proceed to the bottom
      if(((**node).left_node) != NULL) build_tree(&((**node).left_node), node, start, (**node).cpt, th, contrasts, eval_contrast_fun);
      if(((**node).right_node) != NULL) build_tree(&((**node).right_node), node, (**node).cpt+1, end, th, contrasts, eval_contrast_fun);
      
    }
    
  }
 }
  
}

void get_changepoints(cpt_tree_node_t **node, cpts_t *cpts, int start, int end, int min_dist){
  
  if((*node)!= NULL) {
    
    if(((**node).cpt - start + 1 >= min_dist) & (end - (**node).cpt >= min_dist)){
      (*cpts).cpt[(*cpts).n_cpt] = (**node).cpt;
      // NOTE: indice in C starts from 0, so adding 1 to match the indices in R
      (*cpts).index[(*cpts).n_cpt] = ((**node).index)[0]+1;
      (*cpts).n_cpt++;
    }
    
    if((**node).max < (*cpts).min_max) (*cpts).min_max = (**node).max;
    
    get_changepoints(&((**node).left_node), cpts, start, (**node).cpt, min_dist);
    get_changepoints(&((**node).right_node), cpts, (**node).cpt+1, end, min_dist);
    
  }
  
}

int compare_unsigned_int(const void *a, const void *b){
  const int *x = a, *y = b;
  if(*x > *y) return 1;
  else return(*x < *y) ? -1: 0;
}

int compare_cpts_t(const cpts_t *a, const cpts_t *b, int n_obs){
  
  if((*a).n_cpt != (*b).n_cpt) return 1;
  else {
    
    char *tmp = Calloc(n_obs, char);
    memset(tmp, 0,  n_obs * sizeof(char));
    
    int i = 0;
    int are_different = 0;
    
    for(i=0; i<(*a).n_cpt; i++) tmp[(*a).cpt[i]] = 1;
    
    i=0;
    
    while((i<(*a).n_cpt) & (are_different==0)) {
      
      if(tmp[(*b).cpt[i]] != 1)  are_different = 1;
      i++;
       
    } 
    
    Free(tmp);
    
    return are_different;
    
  }
  
}


solution_path_t *solution_path(contrasts_t *contrasts, eval_contrast_fun_t eval_contrast_fun, int min_dist){
  
  //create the solution path
  solution_path_t *solution_path = Calloc(1, solution_path_t);
  (*solution_path).cpts = Calloc(0, cpts_t);
  int len = 0, allocated_len = 0, cpts_not_eqal = 1;
  
  
  cpts_t tmp_cpts;
  double th = 0.0;

  
  tmp_cpts.cpt = Calloc((*contrasts).n_obs, int);
  tmp_cpts.index = Calloc((*contrasts).n_obs, int);
  
  //build the initial tree
  cpt_tree_node_t * root = NULL;
  cpt_tree_node_t * null_node = NULL;

  build_tree(&root, &null_node, 1, (*contrasts).n_obs, th, contrasts, eval_contrast_fun);
 
  while(root != NULL){
    
    //reallocate space for cpts if necessary
    if(allocated_len == len){
      allocated_len += CPTS_LEN_STEP;
      (*solution_path).cpts = Realloc((*solution_path).cpts, allocated_len, cpts_t);
    }
    
    
    // zero the tmp variable
    tmp_cpts.n_cpt = 0;
    tmp_cpts.min_max = DBL_MAX;
    
    // get changepoints (as well as other info, e.g. indices of the subinterval) from the tree
    get_changepoints(&root, &tmp_cpts, 1,  (*contrasts).n_obs, min_dist);
    th = tmp_cpts.min_max;
    
    // check if the new changepoints are the same as the previous ones
    if(len > 0) cpts_not_eqal = compare_cpts_t(&tmp_cpts, &((*solution_path).cpts[len-1]), (*contrasts).n_obs);
  
    if(cpts_not_eqal != 0){
      
      //allocate memory for cpt locations and other relevant info
      (*solution_path).cpts[len].cpt = Calloc(tmp_cpts.n_cpt, int);
      memcpy((*solution_path).cpts[len].cpt, tmp_cpts.cpt, tmp_cpts.n_cpt * sizeof(int));
      (*solution_path).cpts[len].index = Calloc(tmp_cpts.n_cpt, int);
      memcpy((*solution_path).cpts[len].index, tmp_cpts.index, tmp_cpts.n_cpt * sizeof(int));
      (*solution_path).cpts[len].n_cpt = tmp_cpts.n_cpt;
      (*solution_path).cpts[len].min_max = tmp_cpts.min_max;

      len++;
    }

    //reconstruct the tree with a larger threshold
    build_tree(&root, &null_node, 1, (*contrasts).n_obs, th, contrasts, eval_contrast_fun);

 }
  
  (*solution_path).n_th = len;
  //clean_up  
  destroy_tree(&root);
  Free(tmp_cpts.cpt);
  
  return solution_path;
  
}

void destroy_solution_path(solution_path_t **solution_path){
  
  if( (*solution_path) != NULL){
    
    for(int i=0; i<(**solution_path).n_th; i++) Free((**solution_path).cpts[i].cpt);
    Free((**solution_path).cpts);
    Free((**solution_path).th);
    Free(*solution_path);
    
  }
  
  *solution_path = NULL;
  
}

void destroy_tree(cpt_tree_node_t **node){
  
  if( (*node) != NULL){
    
    if( (**node).left_node != NULL) destroy_tree(&((**node).left_node));
    if( (**node).right_node != NULL) destroy_tree(&((**node).right_node));
    if( (**node).index != NULL) Free((**node).index);
    Free(*node);
    
  }
  
  *node = NULL;

}

