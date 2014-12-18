
library(plink2R)

root <- "data"

system.time({
   r1 <- read_plink(root)
})
dim(r1$bed)
sum(is.na(r1$bed))
sum(r1$bed == 3)
max(r1$bed)

system.time({
   r2 <- read_plink(root, impute="random")
})
dim(r2$bed)
sum(is.na(r2$bed))
sum(r2$bed == 3, na.rm=TRUE)
max(r2$bed, na.rm=TRUE)

