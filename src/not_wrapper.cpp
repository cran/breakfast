#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

extern "C" SEXP not_r_wrapper(SEXP x, SEXP intervals, SEXP method, SEXP contrast_type, SEXP parallel, SEXP augmented);

//' @keywords internal
// [[Rcpp::export]]
SEXP call_not_r_wrapper(SEXP x, SEXP intervals, SEXP method, SEXP contrast_type, SEXP parallel, SEXP augmented) {
    return not_r_wrapper(x, intervals, method, contrast_type, parallel, augmented);
}
