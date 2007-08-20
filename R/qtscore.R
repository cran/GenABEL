"qtscore" <-
function(formula,data,snpsubset,idsubset,strata,trait.type="gaussian",times=1,quiet=FALSE,bcast=10,clambda=TRUE,propPs=1.0) {
  	if (class(data)!="gwaa.data") {
		stop("wrong data class: should be gwaa.data")
  	}
	if (!missing(snpsubset)) data <- data[,snpsubset]
	if (!missing(idsubset)) data <- data[idsubset,]
	if (missing(strata)) {nstra=1; strata <- rep(0,data@gtdata@nids)}
	ttargs <- c("gaussian","binomial")
	if (!(match(trait.type,ttargs,nomatch=0)>0)) {
		out <- paste("trait.type argument should be one of",ttargs,"\n")
		stop(out)
	}
	if (trait.type=="gaussian") fam <- gaussian()
	if (trait.type=="binomial") fam <- binomial()

	if (!missing(data)) attach(data@phdata,pos=2,warn.conflicts=FALSE)
	if (class(formula) == "formula") {
		mf <- model.frame(formula,data,na.action=na.omit,drop.unused.levels=TRUE)
		y <- model.response(mf)
#		sdy <- sd(y)
#		meany <- mean(y)
#		if (trait.type=="gaussian") y <- (y-mean(y))/sd(y)
		desmat <- model.matrix(formula,mf)
		lmf <- glm.fit(desmat,y,family=fam)
#		iniest <- lmf$coeff
		mids <- rownames(data@phdata) %in% rownames(mf)
		if (trait.type=="binomial") {
			bin <- 1
		} else {
			bin <- 0
		}
		resid <- lmf$resid
	} else if (class(formula) == "numeric" || class(formula) == "integer" || class(formula) == "double") {
		y <- formula
		mids <- (!is.na(y))
		y <- y[mids]
		resid <- y
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

	} else {
		stop("formula argument must be a formula or one of (numeric, integer, double)")
	}
	if (!missing(data)) detach(data@phdata)
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

	tmeas <- as.logical(mids)
	strata <- strata[tmeas]

	if (any(tmeas == FALSE)) {
		if (!quiet) warning(paste(sum(!tmeas),"people (out of",length(tmeas),") excluded because they have trait or covariate missing\n"),immediate. = TRUE)
		data <- data[tmeas,]
	}

	lenn <- data@gtdata@nsnps;
	out <- list()
	for (j in c(1:(times+1*(times>1)))) {
		if (j>1) resid <- sample(resid,replace=FALSE)
		chi2 <- .C("qtscore",as.raw(data@gtdata@gtps),as.double(resid),as.integer(bin),as.integer(data@gtdata@nids),as.integer(data@gtdata@nsnps), as.integer(nstra), as.integer(strata), chi2 = double(6*data@gtdata@nsnps), PACKAGE="GenABEL")$chi2
		if (any(data@gtdata@chromosome=="X")) {
		  ogX <- data@gtdata[,data@gtdata@chromosome=="X"]
		  sxstra <- strata; sxstra[ogX@male==1] <- strata[ogX@male==1]+nstra
		  chi2X <- .C("qtscore",as.raw(ogX@gtps),as.double(resid),as.integer(bin),as.integer(ogX@nids),as.integer(ogX@nsnps), as.integer(nstra*2), as.integer(sxstra), chi2 = double(6*ogX@nsnps), PACKAGE="GenABEL")$chi2
		  revec <- (data@gtdata@chromosome=="X")
		  revec <- rep(revec,6)
		  chi2 <- replace(chi2,revec,chi2X)
		  rm(ogX,chi2X,revec);gc(verbose=FALSE)
		}
		if (j == 1) {
			chi2.1df <- chi2[1:lenn];
			out$chi2.1df <- chi2.1df
			chi2.2df <- chi2[(lenn+1):(2*lenn)];
			out$chi2.2df <- chi2.2df
			actdf <- chi2[(2*lenn+1):(3*lenn)];
			if (lenn<=10 && !is.numeric(clambda)) {
				lambda <- list()
				lambda$estimate <- NA
				lambda$se <- NA
				chi2.c1df <- chi2.1df;
			} else {
				if (is.numeric(clambda)) {
					lambda <- list()
					lambda$estimate <- clambda
					lambda$se <- NA
					chi2.c1df <- chi2.1df/lambda$estimate;
				} else {
					lambda <- estlambda(chi2.1df,plot=FALSE,prop=propPs)
					def <- 1/lambda$estimate
					if (def > 1 && clambda) {
						chi2.c1df <- chi2.1df;
					} else {
						chi2.c1df <- def*chi2.1df;
					}
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
