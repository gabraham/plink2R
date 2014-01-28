plink2R
=======

plink2R natively reads PLINK (http://pngu.mgh.harvard.edu/~purcell/plink)
BED/BIM/FAM files into R.

Missing genotypes are imputed by assigning missing values the per-SNP average.

Depends on Rcpp, RcppEigen

Example:
   ```
   library(plink2R)
   dat <- read_plink("data")
   dim(dat$bed)
   dim(dat$fam)
   dim(dat$bim)
   ```

Also see the file plink2.R (required glmnet and doMC)

