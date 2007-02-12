"emp.qtscore" <-
function(formula,data,strata,snpsubset,idsubset,times=100,bcast=25) {
  	if (class(data)!="gwaa.data") {
		stop("wrong data class: should be gwaa.data")
  	}
  	if (!is.character(formula)) stop("formula must be character or formula object (apply \"s)")
	if (!missing(snpsubset)) data <- data[,snpsubset]
	if (!missing(idsubset)) data <- data[idsubset,]
	if (missing(strata)) {nstra=1; strata <- rep(0,data@gtdata@nids)}

	avrs <- all.vars(as.formula(formula))
	if (!any(avrs=="CRSNP")) stop("formula must contain CRSNP variable to be replaced with the analysis SNPs")
	if (length(strata)!=data@gtdata@nids) stop("Strata variable and the data do not match in length")
	if (any(is.na(strata))) stop("Strata variable contains NAs")
	if (any(strata!=0)) {
		olev <- levels(as.factor(strata))
		nstra <- length(olev)
		tstr <- strata
		for (i in 0:(nstra-1)) tstr <- replace(tstr,(strata==olev[i+1]),i)
		strata <- tstr
		rm(tstr)
	}
	nstra <- length(levels(as.factor(strata)))

	attach(data@phdata,warn.conflicts=FALSE)
	tmeas <- !is.na(get(avrs[1]))
	avrs <- avrs[avrs!="CRSNP"]
	if (length(avrs)>=2) for (i in 2:length(avrs)) tmeas <- tmeas*(!is.na(get(avrs[i])))
	tmeas <- as.logical(tmeas)
	resid <- get(avrs[1])[tmeas]
	strata <- strata[tmeas]
	if (length(avrs)>=1) for (i in 1:length(avrs)) assign(avrs[i],get(avrs[i])[tmeas])
	detach(data@phdata)

	if (any(tmeas == FALSE)) {
		warning(paste(sum(!tmeas),"people (out of",length(tmeas),") excluded because they have trait or covariate missing\n"),immediate. = TRUE)
		data <- data[tmeas,]
		gc(verbose=FALSE)
	}
	if (length(avrs) > 1) {
		CRSNP <- rep(0,data@gtdata@nids)
		resid <- lm(as.formula(formula),data=data@phdata)$residuals
		rm(CRSNP);gc(verbose=FALSE)
	} 

	a <- .C("qtscore",as.raw(data@gtdata@gtps),as.double(resid),as.integer(data@gtdata@nids),as.integer(data@gtdata@nsnps), as.integer(nstra), as.integer(strata), chi2 = double(6*data@gtdata@nsnps), PACKAGE="GenABEL")$chi2
        lena <- data@gtdata@nsnps
	oP1 <- (1. - pchisq(a[1:lena],1))
	oP2 <- (1. - pchisq(a[(lena+1):(2*lena)],df=(a[(2*lena+1):(3*lena)])))
	effB <- a[(3*lena+1):(lena*4)]
	effAB <- a[(4*lena+1):(lena*5)]
	effBB <- a[(5*lena+1):(lena*6)]
	rm(a);gc(verbose=FALSE)

	nsnps <- data@gtdata@nsnps;
	prob1 <- rep(0,nsnps)
	prob2 <- rep(0,nsnps)
	med1df <- 0
	med2df <- 0
	for (i in 1:times) {
		residx <- sample(resid, replace=FALSE)
		chi2 <- .C("qtscore",as.raw(data@gtdata@gtps),as.double(residx),as.integer(data@gtdata@nids),as.integer(nsnps), as.integer(nstra), as.integer(strata), chi2 = double(6*nsnps), PACKAGE="GenABEL")$chi2
		actdf <- chi2[(2*nsnps+1):(nsnps*3)]
		P1df <- 1. - pchisq(chi2[1:nsnps],1)
		P2df <- 1. - pchisq(chi2[(nsnps+1):(nsnps*2)],df=actdf)
		med1df <- med1df + median(chi2[1:nsnps])
		med2df <- med2df + median(chi2[(nsnps+1):(nsnps*2)])
		mx1 <- min(P1df)
		mx2 <- min(P2df)
		prob1 <- prob1 + 1.*(oP1>mx1)
		prob2 <- prob2 + 1.*(oP2>mx2)
		if (bcast && (i/bcast == round(i/bcast))) print(i)
		gc(verbose=FALSE)
	}
	rm(P1df,P2df,actdf);gc(verbose=FALSE)
	med1df <- med1df/times
	med2df <- med2df/times
        prob1 <- prob1/times
        prob2 <- prob2/times
	if (any(prob1<=0)) {prob1 <- replace(prob1,(prob1<=0),1/(times+1))}
	if (any(prob2<=0)) {prob2 <- replace(prob2,(prob2<=0),1/(times+1))}
	out <- list(P1df = prob1, P2df=prob2, medChi1df = med1df, medChi2df = med2df, name = data@gtdata@snpnames, ids = data@gtdata@idnames, formula = match.call(), family = "qtscore on 1 and 2 d.f.", map = data@gtdata@map, chromosome = data@gtdata@chromosome, effB=effB, effAB=effAB, effBB=effBB)
	rm(prob1,prob2,med1df,med2df);gc(verbose=FALSE)
	class(out) <- "scan.gwaa"
	out
}

