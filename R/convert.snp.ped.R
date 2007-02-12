"convert.snp.ped" <-
function(pedfile,mapfile,outfile,bcast=10000) {
	imap <- read.table(mapfile,header=TRUE)
	iped <- read.table(pedfile,header=FALSE)
	ofile <- file(outfile,"w")
#	ids <- scan(file=ifile,what=character(),nlines=1,quiet=TRUE)
	ids <- as.character(iped[,2])
	nids <- length(ids)
#	mnams <- scan(file=ifile,what=character(),nlines=1,quiet=TRUE)
	mnams <- as.character(imap[,3])
	nmrk <- length(mnams)
	if ((dim(iped)[2]-6)/2 != nmrk) stop("Number of markers in map and ped-files are not equal!")
#	chrom <- scan(file=ifile,what=character(),nlines=1,quiet=TRUE)
	chrom <- as.character(imap[,1])
#	pos <- scan(file=ifile,what=double(),nlines=1,quiet=TRUE)
	pos <- imap[,2]
	nbytes <- ceiling(nids/4)
	cat(file=ofile,ids,"\n")
	cat(file=ofile,mnams,"\n")
	cat(file=ofile,chrom,"\n")
	cat(file=ofile,pos,"\n")
	rdta <- raw(nbytes)
	for (i in 1:nmrk) {
#		gtin <- scan(file=ifile,what=integer(),nlines=1,quiet=TRUE)
		strt <- 6+(i*2-1)
		fin <- strt+1
		a1 <- iped[,strt]
		l1 <- levels(as.factor(a1))
		a2 <- iped[,fin]
		l2 <- levels(as.factor(a2))
		alleles <- sort(as.numeric(unique(c(l1,l2))))
		if (alleles[1] == 0) alleles<-alleles[2:length(alleles)]
		if (length(alleles)<2) alleles[2]=100
		g1 <- alleles[1]*2
		g2 <- alleles[1]+alleles[2]
		g3 <- alleles[2]*2
		gtin <- a1+a2
		gtin <- replace(gtin,(gtin==g1),1)
		gtin <- replace(gtin,(gtin==g2),2)
		gtin <- replace(gtin,(gtin==g3),3)
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
	close(ofile)
}

