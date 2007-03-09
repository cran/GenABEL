"qtscore" <-
function(formula,data,snpsubset,idsubset,strata,trait.type,times=1,quiet=FALSE,bcast=10,clambda=TRUE) {
  	if (class(data)!="gwaa.data") {
		stop("wrong data class: should be gwaa.data")
  	}
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

	attach(data@phdata,pos=2,warn.conflicts=FALSE)
	tmeas <- !is.na(get(avrs[1],pos=2))
	avrs <- avrs[avrs!="CRSNP"]
	if (length(avrs)>=2) for (i in 2:length(avrs)) tmeas <- tmeas*(!is.na(get(avrs[i],pos=2)))
	tmeas <- as.logical(tmeas)
	resid <- get(avrs[1],pos=2)[tmeas]
	strata <- strata[tmeas]
	if (length(avrs)>=1) for (i in 1:length(avrs)) assign(avrs[i],get(avrs[i],pos=2)[tmeas])
	detach(data@phdata)

	if (length(unique(resid))==1) stop("trait is monomorphic")
	if (length(unique(resid))==2) bin <- 1 else bin <- 0
	if (!missing(trait.type)) {
		if (trait.type=="binomial") {
			if (bin == 0) stop("more then two levels in trait")
			bin <- 1
		} else if (trait.type=="gaussian") {
			if (bin == 1) warning("binary traits is analysed as gaussian")
			bin <- 0
		} else {
			stop("trait type should be either binomial or gaussian")
		}
	}
	if (bin == 1) trait.type="binomial" else if (bin == 0) trait.type="gaussian"

	if (any(tmeas == FALSE)) {
		if (!quiet) warning(paste(sum(!tmeas),"people (out of",length(tmeas),") excluded because they have trait or covariate missing\n"),immediate. = TRUE)
		data <- data[tmeas,]
	}
	if (length(avrs) > 1) {
		if (bin) {
			CRSNP <- rep(0,data@gtdata@nids)
			resid <- glm(as.formula(formula),data=data@phdata,family=binomial())$residuals
			resid <- exp(resid)/(1+exp(resid))
		} else {
			CRSNP <- rep(0,data@gtdata@nids)
			resid <- lm(as.formula(formula),data=data@phdata)$residuals
		}
	} 

	lenn <- data@gtdata@nsnps;
	for (j in c(1:(times+1*(times>1)))) {
		if (j>1) resid <- sample(resid,replace=FALSE)
		chi2 <- .C("qtscore",as.raw(data@gtdata@gtps),as.double(resid),as.integer(bin),as.integer(data@gtdata@nids),as.integer(data@gtdata@nsnps), as.integer(nstra), as.integer(strata), chi2 = double(6*data@gtdata@nsnps), PACKAGE="GenABEL")$chi2
		if (any(data@gtdata@chromosome=="X")) {
		  ogX <- data@gtdata[,data@gtdata@chromosome=="X"]
		  chi2X <- .C("qtscore",as.raw(ogX@gtps),as.double(resid),as.integer(bin),as.integer(ogX@nids),as.integer(ogX@nsnps), as.integer(2), as.integer(ogX@male), chi2 = double(6*ogX@nsnps), PACKAGE="GenABEL")$chi2
		  revec <- (data@gtdata@chromosome=="X")
		  revec <- rep(revec,6)
		  chi2 <- replace(chi2,revec,chi2X)
		  rm(ogX,chi2X,revec);gc(verbose=FALSE)
		}
		if (j == 1) {
			chi2.1df <- chi2[1:lenn];
			chi2.2df <- chi2[(lenn+1):(2*lenn)];
			actdf <- chi2[(2*lenn+1):(3*lenn)];
			if (lenn<=10) {
				lambda <- list()
				lambda$estimate <- NA
				lambda$se <- NA
				chi2.c1df <- chi2.1df;
			} else {
				lambda <- estlambda(chi2.1df,plot=FALSE,prop=0.9)
				def <- 1/lambda$estimate
				if (def > 1 && clambda) {
					chi2.c1df <- chi2.1df;
				} else {
					chi2.c1df <- def*chi2.1df;
				}
			}
			effB <- chi2[(3*lenn+1):(lenn*4)]
			effAB <- chi2[(4*lenn+1):(lenn*5)]
			effBB <- chi2[(5*lenn+1):(lenn*6)]
			if (times>1) {
				pr.1df <- rep(0,lenn)
				pr.2df <- rep(0,lenn)
				pr.c1df <- rep(0,lenn)
			}
		} else {
			th1 <- max(chi2[1:lenn])
			pr.1df <- pr.1df + 1*(chi2.1df < th1)
			pr.2df <- pr.2df + 1*(chi2.2df < max(chi2[(lenn+1):(2*lenn)]))
			pr.c1df <- pr.c1df + 1*(chi2.c1df < th1)
			if (!quiet && ((j-1)/bcast == round((j-1)/bcast))) {
				cat("\b\b\b\b\b\b",round((100*(j-1)/times),digits=2),"%",sep="")
				flush.console()
			}
		}
	}
	if (times > bcast) cat("\n")

	out <- list()
	if (times>1) {
		out$P1df <- pr.1df/times
		out$P1df <- replace(out$P1df,(out$P1df==0),1/(1+times))
		out$P2df <- pr.2df/times
		out$P2df <- replace(out$P2df,(out$P2df==0),1/(1+times))
		out$Pc1df <- pr.c1df/times
		out$Pc1df <- replace(out$Pc1df,(out$Pc1df==0),1/(1+times))
	} else {
		out$P1df <- 1. - pchisq(chi2.1df,1)
		out$P2df <- 1. - pchisq(chi2.2df,actdf)
		out$Pc1df <- 1. - pchisq(chi2.c1df,1)
	}
	out$lambda <- lambda
	out$effB <- effB
	out$effAB <- effAB
	out$effBB <- effBB
	out$snpnames <- data@gtdata@snpnames
	out$map <- data@gtdata@map
	out$chromosome <- data@gtdata@chromosome
	out$idnames <- data@gtdata@idnames
	out$formula <- match.call()
	out$family <- paste("score test for association with trait type",trait.type)
	class(out) <- "scan.gwaa"
	out
}
