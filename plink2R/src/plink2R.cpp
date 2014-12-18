
#include "plink2R.h"
#include "data.hpp"

// [[Rcpp::depends(BH)]]
SEXP read_plink(SEXP _bedfile, SEXP _famfile,
   SEXP _impute, SEXP _verbose) 
{
   std::string bedfile = as<std::string>(_bedfile);
   std::string famfile = as<std::string>(_famfile);
   int impute = as<int>(_impute);
   bool verbose = as<bool>(_verbose);

   if(verbose)
      std::cout << "[read_plink] bedfile: " << bedfile <<
	 " famfile:" << famfile << std::endl;
   Data data(bedfile.c_str(), famfile.c_str(), verbose);
   data.read_bed(impute);

   if(impute)
   {
      return wrap(data.X);
   }
   else
   {
      NumericMatrix X(wrap(data.X));
      NumericVector X2 = wrap(ifelse(X == PLINK_NA, NA_REAL, X));
      X2.attr("dim") = Dimension(X.nrow(), X.ncol());
      return X2;
   }
}

