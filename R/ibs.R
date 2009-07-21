"ibs" <- 
function (data,snpsubset,idsubset,cross.idsubset,weight="no",snpfreq) {
# idsubset, cross.idsubset: should be real names, not indexes!
	if (is(data,"gwaa.data")) data <- data@gtdata
	if (!is(data,"snp.data")) stop("The data argument must be of snp.data-class or gwaa.data-class")
	if (!missing(snpsubset)) data <- data[,snpsubset]
	if (missing(idsubset) && !missing(cross.idsubset)) stop("cross.idsubset arg cannot be used (idsubset missing)",immediate. = TRUE)
	if (!missing(snpfreq)) {
		if (length(snpfreq) != data@nsnps) stop("snpfreq argument not equal in length to the number of SNPs in data")
		if (any(snpfreq<0.) || any(snpfreq>1.)) stop("snpfreq argument: frequencies out of [0,1]")
		if (!is(snpfreq,"numeric")) stop("snpfreq argument: non-numeric class")
	} else {
		snpfreq <- summary(data)[,"Q.2"]
	}
	if (weight=="no") {
		homodiag <- hom(data)[,"Hom"]
	} else {
		homodiag <- 0.5+hom(data,snpfreq=snpfreq)[,"F"]
	}
	ibs.C.option <- 0
	if (!missing(idsubset) && !(is.numeric(idsubset) || is.logical(idsubset) || is.character(idsubset))) stop("idsubset must be numeric, logical, or character")
	if (!missing(cross.idsubset) && !(is.numeric(cross.idsubset) || is.logical(cross.idsubset) || is.character(cross.idsubset))) stop("cross.idsubset must be numeric, logical, or character")
	if (!missing(idsubset) && (is.numeric(idsubset) || is.logical(idsubset))) idsubset <- data@idnames[idsubset]
	if (!missing(cross.idsubset) && (is.numeric(cross.idsubset) || is.logical(cross.idsubset))) cross.idsubset <- data@idnames[cross.idsubset]
	if (!missing(idsubset) && !missing(cross.idsubset)) {
		idset1 <- idsubset
		idset2 <- cross.idsubset
		if (any(idset1 %in% idset2)) stop("idsubset and cross.idsubset should not overlap!")
		idsorder <- c(idset1,idset2)
		homodiag <- homodiag[match(idsorder,data@idnames)]
		if (length(idsorder) != data@nids) data <- data[idsorder,]
		if (any(idsorder!=data@idnames)) data <- data[idsorder,]
		ibs.C.option <- 1
	} else if (!missing(idsubset) & missing(cross.idsubset)) {
		idset1 <- idsubset
		idset2 <- idsubset
		idsorder <- idsubset
		homodiag <- homodiag[match(idsorder,data@idnames)]
		data <- data[idsorder,]
	} else if (missing(idsubset) & missing(cross.idsubset)) {
		idset1 <- data@idnames
		idset2 <- data@idnames
	} else {
		stop("can not be: impossible combination of idsubset and cross.idsubset")
	}
	gc()
	idset1.num <- which(data@idnames %in% idset1)
	idset2.num <- which(data@idnames %in% idset2)
	wargs <- c("no","freq")
	if (!(match(weight,wargs,nomatch=0)>0)) {
		out <- paste("weight argument should be one of",wargs,"\n")
		stop(out)
	}
#	if (npairs > 2500000) stop("Too many pairs to evaluate... Stopped")
	if (weight == "no") option = 0
	if (weight == "freq") option = 1
	if (ibs.C.option==1) {
		sout <- .C("ibspar", as.raw(data@gtps), as.integer(data@nids), as.integer(data@nsnps), as.integer(length(idset1.num)), as.integer(idset1.num-1), as.integer(length(idset2.num)), as.integer(idset2.num-1), as.double(snpfreq), as.integer(option), sout = double(2*length(idset1.num)*length(idset2.num)), PACKAGE="GenABEL")$sout
		out <- list()
		out$ibs <- sout[1:(length(idset1.num)*length(idset2.num))]
		out$num <- sout[(length(idset1.num)*length(idset2.num)+1):(length(idset1.num)*length(idset2.num)*2)]
		dim(out$ibs) <- c(length(idset2.num),length(idset1.num))
		dim(out$num) <- c(length(idset1.num),length(idset2.num))
		rownames(out$ibs) <- idset2
		colnames(out$ibs) <- idset1
		rownames(out$num) <- idset1
		colnames(out$num) <- idset2
	} else if (ibs.C.option==0) {
		out <- .C("ibsnew", as.raw(data@gtps), as.integer(data@nids), as.integer(data@nsnps), as.double(snpfreq), as.integer(option), sout = double(data@nids*data@nids), PACKAGE="GenABEL")$sout
		dim(out) <- c(length(idset2.num),length(idset1.num))
		out <- t(out)
		diag(out) <- homodiag
		rownames(out) <- idset1
		colnames(out) <- idset2
	} else {
		stop("can not be: incorrect ibs.C.option")
	}
	out
}
