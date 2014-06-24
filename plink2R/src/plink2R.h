#ifndef _plink2R_RCPP_READ_PLINK
#define _plink2R_RCPP_READ_PLINK

#include <RcppEigen.h>

#define IMPUTE_NONE 0
#define IMPUTE_AVG 1
#define IMPUTE_RANDOM 2

using namespace Rcpp;
using namespace Eigen;

RcppExport SEXP read_plink(SEXP bedfile, SEXP famfile, SEXP impute,
   SEXP verbose);

#endif
