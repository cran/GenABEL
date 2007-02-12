"convert.snp.text" <-
function(infile,outfile,bcast=10000) {
	ifile <- file(infile,"r")
	ofile <- file(outfile,"w")
	ids <- scan(file=ifile,what=character(),nlines=1,quiet=TRUE)
	nids <- length(ids)
	mnams <- scan(file=ifile,what=character(),nlines=1,quiet=TRUE)
	nmrk <- length(mnams)
	chrom <- scan(file=ifile,what=character(),nlines=1,quiet=TRUE)
	pos <- scan(file=ifile,what=double(),nlines=1,quiet=TRUE)
	nbytes <- ceiling(nids/4)
	cat(file=ofile,ids,"\n")
	cat(file=ofile,mnams,"\n")
	cat(file=ofile,chrom,"\n")
	cat(file=ofile,pos,"\n")
	rdta <- raw(nbytes)
	for (i in 1:nmrk) {
		gtin <- scan(file=ifile,what=integer(),nlines=1,quiet=TRUE)
		gchk <- (gtin==0 | gtin==1 | gtin ==2 | gtin ==3)
		if (!all(gchk)) {
			cat("Wrong genotype codes:\nCODE\tID\tSNP\n")
			wlst <- which(!gchk)
			for (j in 1:length(wlst)) cat(gtin[wlst[j]],"\t",ids[wlst[j]],"\t",mnams[i],"\n")
			stop("execution terminated")
		}
		rdta <- put.snps(gtin)
		cat(file=ofile,rdta,"\n")
		if (bcast && round(i/bcast) == i/bcast) cat("Converted",i,"records...\n")
	}
	close(ifile)
	close(ofile)
}

