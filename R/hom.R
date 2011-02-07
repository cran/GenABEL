"hom" <- 
		function(data,snpsubset,idsubset,snpfreq,n.snpfreq=1000) {
	if (is(data,"gwaa.data")) {
		data <- data@gtdata
	}
	if (!is(data,"snp.data")) {
		stop("data argument should have gwaa.data-class or snp.data-class")
	}
	
	if (!missing(snpsubset)) data <- data[,snpsubset]
	useFreq <- 0
	if (!missing(snpfreq)) {
		if (length(snpfreq) != data@nsnps) stop("snpfreq argument not equal in length to the number of SNPs in data")
		if (any(snpfreq<0) || any(snpfreq>1.)) stop("snpfreq argument: frequencies out of [0,1]")
		if (!is(snpfreq,"numeric")) stop("snpfreq argument: non-numeric class")
		useFreq <- 1
	} else {
		snpfreq <- rep(-1.,data@nsnps)
	}
	if (!missing(idsubset)) data <- data[idsubset,]
	
	totsnps <- data@nsnps
	totids <- data@nids
	if (!is(n.snpfreq,"numeric") || min(n.snpfreq)<0) stop("n.snpfreq should be positive numeric")
	if (length(n.snpfreq)==1) {
		n.snpfreq <- rep(n.snpfreq,data@nsnps)
	} else {
		if (length(n.snpfreq) != data@nsnps) stop("length mismatch in n.snpfreq")
	}
	out <- .C("hom",as.raw(data@gtps),as.integer(data@nids),as.integer(data@nsnps),as.double(snpfreq),as.double(n.snpfreq),as.integer(useFreq),sout = double(5*data@nids), PACKAGE="GenABEL")$sout
	dim(out) <- c(data@nids,5)
	F <- (out[,3]-out[,4])/(out[,1]-out[,4])
	out <- cbind(out,F)
	out[,3] <- out[,3]/out[,1]
	out[,4] <- out[,4]/out[,1]
	out[,5] <- out[,5]/out[,2]
	colnames(out) <- c("NoMeasured","NoPoly","Hom","E(Hom)","Var","F")
	out <- as.data.frame(out,stringsAsFactors=FALSE)
	rownames(out) <- data@idnames
	out
}
