"qtscore" <-
function(formula,data,strata,snpsubset,idsubset,quiet=FALSE) {
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
		if (!quiet) warning(paste(sum(!tmeas),"people (out of",length(tmeas),") excluded because they have trait or covariate missing\n"),immediate. = TRUE)
		data <- data[tmeas,]
	}
	if (length(avrs) > 1) {
		CRSNP <- rep(0,data@gtdata@nids)
		resid <- lm(as.formula(formula),data=data@phdata)$residuals
	} 

	chi2 <- .C("qtscore",as.raw(data@gtdata@gtps),as.double(resid),as.integer(data@gtdata@nids),as.integer(data@gtdata@nsnps), as.integer(nstra), as.integer(strata), chi2 = double(6*data@gtdata@nsnps), PACKAGE="GenABEL")$chi2
	if (any(data@gtdata@chromosome=="X")) {
	  ogX <- data@gtdata[,data@gtdata@chromosome=="X"]
	  chi2X <- .C("qtscore",as.raw(ogX@gtps),as.double(resid),as.integer(ogX@nids),as.integer(ogX@nsnps), as.integer(2), as.integer(ogX@male), chi2 = double(6*ogX@nsnps), PACKAGE="GenABEL")$chi2
	  revec <- (data@gtdata@chromosome=="X")
	  revec <- rep(revec,6)
	  chi2 <- replace(chi2,revec,chi2X)
	  rm(ogX,chi2X,revec);gc(verbose=FALSE)
	}

	lenn <- data@gtdata@nsnps;
	chromosome <- data@gtdata@chromosome
	out <- list(P1df=(1. - pchisq(chi2[1:lenn],1)), P2df=(1. - pchisq(chi2[(lenn+1):(lenn*2)],df=chi2[(2*lenn+1):(lenn*3)])), medChi1df = median(chi2[1:lenn]), medChi2df = median(chi2[(lenn+1):(lenn*2)]), name = data@gtdata@snpnames, formula = formula, family = "qtscore with 1 and 2 d.f.", map = data@gtdata@map, chromosome=data@gtdata@chromosome, ids = data@gtdata@idnames, effB=chi2[(3*lenn+1):(lenn*4)], effAB=chi2[(4*lenn+1):(lenn*5)], effBB=chi2[(5*lenn+1):(lenn*6)])
	class(out) <- "scan.gwaa"
	out
}
