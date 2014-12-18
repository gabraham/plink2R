
read_plink <- function(root, snps=NULL, impute=c("none", "avg", "random"), verbose=FALSE)
{
   proot <- path.expand(root)
   impute <- match.arg(impute)
   impute_int <- switch(impute,
      "none"=0L,
      "avg"=1L,
      "random"=2L
   )

   bedfile <- paste(proot, ".bed", sep="")
   famfile <- paste(proot, ".fam", sep="")
   bimfile <- paste(proot, ".bim", sep="")

   bim <- read.table(bimfile, header=FALSE, sep="", stringsAsFactors=FALSE)
   fam <- read.table(famfile, header=FALSE, sep="", stringsAsFactors=FALSE)
   geno <- .Call("read_plink", PACKAGE="plink2R", bedfile, famfile,
      impute_int, verbose)
   colnames(geno) <- bim[,2]
   rownames(geno) <- paste(fam[,1], fam[, 2], sep=":")

   list(bed=geno, fam=fam, bim=bim)
}

