"ccfast" <-
function(y,data,snpsubset,idsubset,quiet=FALSE) {
	if (class(data) != "gwaa.data") stop("wrong type of data argument, must be gwaa.data")
	if (!missing(snpsubset)) data <- data[,snpsubset]
	if (!missing(idsubset)) data <- data[idsubset,]

	attach(data@phdata,warn.conflicts=FALSE)
	cc <- get(y)
	detach(data@phdata)

        if (length(levels(as.factor(cc)))<2) stop("cc status is monomorphic!") 
        if (length(levels(as.factor(cc)))>2) stop("cc status has more then 2 levels!") 
        if (levels(as.factor(cc))[1] != 0 || levels(as.factor(cc))[2] != 1) stop ("cc is case-control status, with 0 as control and 1 as cases. No 0 and/or 1 found in the data")

	if (any(is.na(cc))) {
	  if (!quiet) warning(paste(sum(is.na(cc)),"people (out of",length(cc),") excluded as not having cc status\n"),immediate. = TRUE)
	  vec <- !is.na(cc)
	  data <- data[vec,]
	  cc1 <- cc[!is.na(cc)]
	} else {
	  cc1 <- cc
	}
	rm(cc)

	a <- fcc(data@gtdata,cc1)
        lena <- data@gtdata@nsnps

	out <- list()
	out$P1df <- (1. - pchisq(a[1:lena],1))
	out$medChi1df <- median(a[1:lena])
	out$P2df <- (1. - pchisq(a[(lena+1):(2*lena)],df=(a[(2*lena+1):(3*lena)])))
	out$medChi2df <- median(a[(lena+1):(lena*2)])
	out$effB <- a[(3*lena+1):(lena*4)]
	out$effAB <- a[(4*lena+1):(lena*5)]
	out$effBB <- a[(5*lena+1):(lena*6)]
	rm(a);gc(verbose=FALSE)
	out$name <- data@gtdata@snpnames
	out$formula <- match.call()
	out$family <- "chi-square 1 and 2 d.f."
	out$map <- data@gtdata@map
	out$chromosome <- data@gtdata@chromosome
	out$ids <- data@gtdata@idnames
	class(out) <- "scan.gwaa"
	out
}

