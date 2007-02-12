"emp.ccfast" <-
function(y,data,snpsubset,idsubset,times=100,bcast=25) {
	if (class(data) != "gwaa.data") stop("wrong type of argument, must be gwaa.data")
	if (!missing(snpsubset)) data <- data[,snpsubset]
	if (!missing(idsubset)) data <- data[idsubset,]


	attach(data@phdata,warn.conflicts=FALSE)
	cc <- get(y)
	detach(data@phdata)

        if (length(levels(as.factor(cc)))<2) stop("cc status is monomorphic!") 
        if (length(levels(as.factor(cc)))>2) stop("cc status has more then 2 levels!") 
        if (levels(as.factor(cc))[1] != 0 || levels(as.factor(cc))[2] != 1) stop ("cc is case-control status, with 0 as control and 1 as cases. No 0 and/or 1 found in the data")

	if (any(is.na(cc))) {
	  warning(paste(sum(is.na(cc)),"people (out of",length(cc),") excluded as not having cc status\n"),immediate. = TRUE)
	  data <- data[!is.na(cc),]
	  cc1 <- cc[!is.na(cc)]
	} else {
	  cc1 <- cc
	}
	rm(cc)

	a <- .C("fastcc",as.raw(data@gtdata@gtps),as.integer(cc1),as.integer(data@gtdata@nids),as.integer(data@gtdata@nsnps), chi2 = double(6*data@gtdata@nsnps), PACKAGE="GenABEL")$chi2
        lena <- data@gtdata@nsnps
	oP1 <- (1. - pchisq(a[1:lena],1))
	oP2 <- (1. - pchisq(a[(lena+1):(2*lena)],df=(a[(2*lena+1):(3*lena)])))
	effB <- a[(3*lena+1):(lena*4)]
	effAB <- a[(4*lena+1):(lena*5)]
	effBB <- a[(5*lena+1):(lena*6)]
	rm(a);gc(verbose=FALSE)

	nsnps <- data@gtdata@nsnps
	med1df <- 0
	med2df <- 0
	prob1 <- rep(0,nsnps)
	prob2 <- rep(0,nsnps)
	P1df <- rep(NA,nsnps)
	P2df <- rep(NA,nsnps)
	med1df <- rep(NA,nsnps)
	med2df <- rep(NA,nsnps)

	for (i in 1:times) {
		ccx <- sample(cc1,replace=FALSE)
		a <- .C("fastcc",as.raw(data@gtdata@gtps),as.integer(ccx),as.integer(data@gtdata@nids),as.integer(nsnps), chi2 = double(6*nsnps), PACKAGE="GenABEL")$chi2
		P1df <- (1. - pchisq(a[1:nsnps],1))
		P2df <- (1. - pchisq(a[(nsnps+1):(2*nsnps)],df=(a[(2*nsnps+1):(3*nsnps)])))
		med1df <- med1df + median(a[1:nsnps])
		med2df <- med2df + median(a[(nsnps+1):(nsnps*2)])
		mx1 <- min(P1df)
		mx2 <- min(P2df)
		prob1 <- prob1 + 1.*(oP1>mx1)
		prob2 <- prob2 + 1.*(oP2>mx2)
		if (bcast && (i/bcast == round(i/bcast))) print(i)
		gc(verbose=FALSE)
	}
	rm(ccx,a,P1df,P2df);gc(verbose=FALSE)

	med1df <- med1df/times
	med2df <- med2df/times
        prob1 <- prob1/times
        prob2 <- prob2/times

	if (any(prob1<=0)) {prob1 <- replace(prob1,(prob1<=0),1/(times+1))}
	if (any(prob2<=0)) {prob2 <- replace(prob2,(prob2<=0),1/(times+1))}

	out <- list(P1df = prob1, P2df=prob2, medChi1df = med1df, medChi2df = med2df, name = data@gtdata@snpnames, ids = data@gtdata@idnames, formula = match.call(), family = "chi-square with 1 and 2 d.f.", map = data@gtdata@map, chromosome = data@gtdata@chromosome, effB=effB, effAB=effAB, effBB=effBB)
#	rm(prob1,prob2,med1df,med2df,effB,effAB,effBB);gc(verbose=FALSE)
	class(out) <- "scan.gwaa"
	out
}

