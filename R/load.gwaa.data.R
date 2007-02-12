"load.gwaa.data" <-
function(phenofile = "pheno.dat", genofile = "geno.raw",force = FALSE, makemap=FALSE) {
# check that ID and SEX are correct
	dta <- read.table(phenofile,header=TRUE,as.is=TRUE)
	coln <- names(dta)
	if (!any(names(dta)=="id",na.rm=TRUE)) 
		stop("the filed named \"id\", containing the identifier presented in both pheno- and geno- files was not found in the phenofile")
	if (!any(names(dta)=="sex",na.rm=TRUE)) 
		stop("the column named \"sex\", containing the male identifier was not found in the phenofile")
	a <- table(dta$sex,exclude=NULL)
	if (length(a) > 2)
		stop("column named \"sex\" contains more then 2 codes")
	if (length(a) == 1 && !(names(a)[1] == 0 || names(a)[1] == 1))
		stop("the column named \"sex\" contains 1 code which is neither 0 (=female) or 1 (=male)")
	if (length(a) == 2 && names(a)[1] != 0 && names(a)[2] != 1)
		stop("the column named \"sex\" is not coded as 0=female and 1=male")
	rm(a);gc(verbose=FALSE)
	if (any(table(dta$id)>1))
		stop("there are duplicated IDs in the phenotypic data file")
	if (any(is.na(dta$id)))
		stop("there are missing IDs in the phenotypic data file")
	if (any(is.na(dta$sex)))
		stop("there are missing sexes in the phenotypic data file")
# read in genotypic data
	ifile <- file(genofile,"r")
	ids <- scan(file=ifile,what=character(),nlines=1,quiet=TRUE)
	nids <- length(ids)
	print("ids loaded...")
	mnams <- scan(file=ifile,what=character(),nlines=1,quiet=TRUE)
	print("marker names loaded...")
	chrom <- scan(file=ifile,what=character(),nlines=1,quiet=TRUE)
	chrom <- as.factor(chrom);gc(verbose=FALSE)
	print("chromosome data loaded...")
	pos <- scan(file=ifile,what=double(),nlines=1,quiet=TRUE)
	print("map data loaded...")
	nsnps <- length(mnams)
	nbytes <- ceiling(nids/4)
	rdta <- scan(file=ifile,what=raw(),quiet=TRUE)
	print("genotype data loaded...")
	close(ifile)
	dim(rdta) <- c(nbytes,nsnps)
	rdta <- new("snp.mx",rdta);gc(verbose=FALSE)

# check errors of match between pheno and geno files
	mlst0 <- match(as.character(dta$id),ids)
	for (i in 1:length(mlst0)) {
		cid <- mlst0[i];
		if (is.na(cid)) 
			cat("person with id =",as.character(dta$id)[i],"was not found in genotypic file; excluded\n")
	}
	mlst <- match(ids,as.character(dta$id))
	oerr <- 0
	for (i in 1:length(mlst)) {
		cid <- mlst[i];
		if (is.na(cid)) {
			cat("person with id =",ids[i],"was not found in phenotypic file!!! - FATAL\n")
			oerr <- oerr + 1
		}
	}
	if (oerr) stop("fatal error. update pheno-file")
	newdta <- data.frame(dta[mlst,])
	rm(dta);gc(verbose=FALSE)

	a <- snp.data(nids=nids,rawdata=rdta,idnames=ids,snpnames=mnams,chromosome=chrom,map=pos,male=newdta$sex)
	print("snp.data object created...")
	rm(rdta,ids,mnams,chrom,pos);gc(verbose=FALSE)

#check X chromosome markers
 	if (any(a@chromosome == "X") && any(a@male == 1) && !force) {
		xmrk <- (a@chromosome == "X")
		mlst <- (a@male == 1)
		rxm <- a[mlst,xmrk]
		rm(xmrk,mlst);gc(verbose=FALSE)
		xerr <- 0
		for (snpX in (1:rxm@nsnps)) {
		 xmdta <- as.double(rxm[,snpX])
		 xm <- (xmdta == 1)
		 xm <- replace(xm,is.na(xm),FALSE)
		 if (any(xm,na.rm=TRUE)) {
			if (!xerr) {
				xerr=1
				tmpoe <- matrix(rep(NA,2*sum(xm)),ncol=2)
				tmpoe[,1] <- rxm@idnames[xm]
				tmpoe[,2] <- rep(rxm@snpnames[snpX],sum(xm))
				outerr <- tmpoe
			}
			tmpoe <- matrix(rep(NA,2*sum(xm)),ncol=2)
			tmpoe[,1] <- rxm@idnames[xm]
			tmpoe[,2] <- rep(rxm@snpnames[snpX],sum(xm))
			outerr <- rbind(outerr,tmpoe)
		  }
		}
		rm(rxm);gc(verbose=FALSE)
		if (xerr) {
			colnames(outerr) <- c("ID","SNP")
			cat("Wrong male X genotypes (heterozygous) found in",dim(outerr)[1],"cases\n")
			cat("Error table is saved as the output object\n")
			return(outerr)
		}
	}
	if (force) cat("assignment of gwaa.data object FORCED; X-errors were not checked!\n")

# make map
	if (makemap) {
		cat("increase in map order FORCED\n")
		gsize <- max(a@map[a@chromosome=="1"])/5
		for (i in 2:22) {
			inc <- max(a@map[a@chromosome==as.character(i-1)]) + gsize
			a@map[a@chromosome==as.character(i)] <- a@map[a@chromosome==as.character(i)] + inc
		}
		inc <- max(a@map[a@chromosome=="22"]) + gsize
		a@map[a@chromosome=="X"] <- a@map[a@chromosome=="X"] + inc
	}	

	out <- new("gwaa.data",phdata=newdta,gtdata=a)
	rm(a,newdta);gc(verbose=FALSE)
	out
}

