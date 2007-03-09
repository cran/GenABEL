"Xcheck" <-
function(data) {
	if (class(data) != "snp.data") stop("data argument should be of snpdata-class")
	if (any(data@chromosome != "X")) stop("All markers should be X-linked")
	xerr <- 0
	for (snpX in (1:data@nsnps)) {
	 xmdta <- as.double(data[,snpX])
	 xm <- (xmdta == 1)
	 xm <- replace(xm,is.na(xm),FALSE)
	 if (any(xm,na.rm=TRUE)) {
		if (!xerr) {
			xerr=1
			tmpoe <- matrix(rep(NA,2*sum(xm)),ncol=2)
			tmpoe[,1] <- data@idnames[xm]
			tmpoe[,2] <- rep(data@snpnames[snpX],sum(xm))
			outerr <- tmpoe
		}
		tmpoe <- matrix(rep(NA,2*sum(xm)),ncol=2)
		tmpoe[,1] <- data@idnames[xm]
		tmpoe[,2] <- rep(data@snpnames[snpX],sum(xm))
		outerr <- rbind(outerr,tmpoe)
	  }
	}
	out <- list();
	out$xerr <- xerr
	if (xerr) {
		colnames(outerr) <- c("ID","SNP")
		out$tab <- outerr
	} else {
		out$tab <- NULL
	}
	out
}
