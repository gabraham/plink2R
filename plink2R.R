
library(plink2R)
library(glmnet)
library(doMC)

registerDoMC(cores=10)

root <- "data"

system.time({
   r1 <- read_plink(root)
})
dim(r1$bed)
sum(is.na(r1$bed))
sum(r1$bed == 3)
max(r1$bed)

system.time({
   r2 <- read_plink(root, impute=FALSE)
})
dim(r2$bed)
sum(is.na(r2$bed))
sum(r2$bed == 3, na.rm=TRUE)
max(r2$bed, na.rm=TRUE)

g <- cv.glmnet(r1$bed, r1$fam[, 6] - 1, 
   family="binomial", nfolds=10, type="auc", dfmax=100,
   parallel=TRUE)

max(g$cvm)

