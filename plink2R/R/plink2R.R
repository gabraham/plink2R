
read_plink <- function(root, snps=NULL, impute=TRUE, verbose=FALSE)
{
   proot <- path.expand(root)
   bedfile <- paste(proot, ".bed", sep="")
   famfile <- paste(proot, ".fam", sep="")
   bimfile <- paste(proot, ".bim", sep="")

   bim <- read.table(bimfile, header=FALSE, sep="", stringsAsFactors=FALSE)
   fam <- read.table(famfile, header=FALSE, sep="", stringsAsFactors=FALSE)
   geno <- .Call("read_plink", PACKAGE="plink2R", bedfile, famfile, impute,
	 verbose)
   colnames(geno) <- bim[,2]
   rownames(geno) <- paste(fam[,1], fam[, 2], sep=":")

   list(bed=geno, fam=fam, bim=bim)
}

