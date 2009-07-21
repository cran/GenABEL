"qtscore" <-
function(formula,data,snpsubset,idsubset,strata,trait.type="gaussian",times=1,quiet=FALSE,bcast=10,clambda=TRUE,propPs=1.0,details=TRUE) {
  	if (class(data)!="gwaa.data") {
		stop("wrong data class: should be gwaa.data")
  	}
	checkphengen(data)
	if (!missing(snpsubset)) data <- data[,snpsubset]
	if (!missing(idsubset)) data <- data[idsubset,]
	if (missing(strata)) {nstra=1; strata <- rep(0,data@gtdata@nids)}
	ttargs <- c("gaussian","binomial","guess")
	if (!(match(trait.type,ttargs,nomatch=0)>0)) {
		out <- paste("trait.type argument should be one of",ttargs,"\n")
		stop(out)
	}
	if (trait.type=="guess") {
		if (!missing(data)) attach(data@phdata,pos=2,warn.conflicts=FALSE)
		if (class(formula) == "formula") {
			mf <- model.frame(formula,data,na.action=na.omit,drop.unused.levels=TRUE)
			y <- model.response(mf)
			if (isbinomial(y)) trait.type <- "binomial" else trait.type <- "gaussian"
		} else if (class(formula) == "numeric" || class(formula) == "integer" || class(formula) == "double") {
			y <- formula
			if (isbinomial(y)) trait.type <- "binomial" else trait.type <- "gaussian"
		} else {
			stop("formula argument must be a formula or one of (numeric, integer, double)")
		}
		warning(paste("trait type is guessed as",trait.type))
		if (!missing(data)) detach(data@phdata)
	}
	if (trait.type=="gaussian") fam <- gaussian()
	if (trait.type=="binomial") fam <- binomial()

	if (!missing(data)) attach(data@phdata,pos=2,warn.conflicts=FALSE)
	if (class(formula) == "formula") {
		mf <- model.frame(formula,data,na.action=na.omit,drop.unused.levels=TRUE)
		y <- model.response(mf)
		test.type(y,trait.type)
		desmat <- model.matrix(formula,mf)
		lmf <- glm.fit(desmat,y,family=fam)
		mids <- rownames(data@phdata) %in% rownames(mf)
		if (trait.type=="binomial") {
			resid <- residuals(lmf,type="response")
			resid <- (resid-min(resid))/(max(resid)-min(resid))
		} else {
			resid <- lmf$resid
		}
	} else if (class(formula) == "numeric" || class(formula) == "integer" || class(formula) == "double") {
		y <- formula
		test.type(y,trait.type)
		mids <- (!is.na(y))
		y <- y[mids]
#		if (trait.type=="binomial") {
#			y <- glm(y~rep(1,length(y)),family=binomial)$resid
#		}
		resid <- y
	} else {
		stop("formula argument should be a formula or a numeric vector")
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
		if (!quiet) warning(paste(sum(!tmeas),"observations deleted due to missingness"))
		data <- data[tmeas,]
	}

#	if (trait.type=="binomial" & class(formula) != "formula") bin<-1 else bin <- 0
	if (trait.type=="binomial") bin<-1 else bin<-0
	lenn <- data@gtdata@nsnps;
	out <- list()
	for (j in c(1:(times+1*(times>1)))) {
		if (j>1) resid <- sample(resid,replace=FALSE)
		chi2 <- .C("qtscore_glob",as.raw(data@gtdata@gtps),as.double(resid),as.integer(bin),as.integer(data@gtdata@nids),as.integer(data@gtdata@nsnps), as.integer(nstra), as.integer(strata), chi2 = double(10*data@gtdata@nsnps), PACKAGE="GenABEL")$chi2
		if (any(data@gtdata@chromosome=="X")) {
		  ogX <- data@gtdata[,data@gtdata@chromosome=="X"]
		  sxstra <- strata; sxstra[ogX@male==1] <- strata[ogX@male==1]+nstra
		  chi2X <- .C("qtscore_glob",as.raw(ogX@gtps),as.double(resid),as.integer(bin),as.integer(ogX@nids),as.integer(ogX@nsnps), as.integer(nstra*2), as.integer(sxstra), chi2 = double(10*ogX@nsnps), PACKAGE="GenABEL")$chi2
		  revec <- (data@gtdata@chromosome=="X")
		  revec <- rep(revec,6)
		  chi2 <- replace(chi2,revec,chi2X)
		  rm(ogX,chi2X,revec);gc(verbose=FALSE)
		}
		if (j == 1) {
			chi2.1df <- chi2[1:lenn];
			chi2.1df[abs(chi2.1df+999.99)<1.e-8] <- 0 #NA
			out$chi2.1df <- chi2.1df
			chi2.2df <- chi2[(lenn+1):(2*lenn)];
			chi2.2df[abs(chi2.2df+999.99)<1.e-8] <- 0 #NA
			actdf <- chi2[(2*lenn+1):(3*lenn)];
			actdf[abs(actdf+999.99)<1.e-8] <- 1.e-16 #NA
#			out$actdf <- actdf
			out$chi2.2df <- chi2.2df
			z0 <- chi2[(7*lenn+1):(8*lenn)];
			z0[abs(z0+999.99)<1.e-8] <- 0 #NA
#			out$z0 <- z0
			z2 <- chi2[(8*lenn+1):(9*lenn)];
			z2[abs(z2+999.99)<1.e-8] <- 0 #NA
#			out$z2 <- z2
			rho <- chi2[(9*lenn+1):(10*lenn)];
			rho[abs(rho+999.99)<1.e-8] <- 0 #NA
#			rho <- abs(rho)
#			out$rho <- rho

			lambda <- list()
			if (is.logical(clambda)) {
				if (lenn<10) {
					warning("no. observations < 10; Lambda set to 1")
					lambda$estimate <- 1.0
					lambda$se <- NA
				} else {
					if (lenn<100) warning("Number of observations < 100, Lambda estimate is unreliable")
					lambda <- estlambda(chi2.1df,plot=FALSE,prop=propPs)
					if (lambda$estimate<1.0 && clambda==TRUE) {
						warning("Lambda estimated < 1, set to 1")
						lambda$estimate <- 1.0
						lambda$se <- NA
					}
				}
			} else {
				if (is.numeric(clambda)) {
					lambda$estimate <- clambda
					lambda$se <- NA
				} else if (is.list(clambda)) {
					if (any(is.na(match(c("estimate","se"),names(clambda)))))
						stop("when clambda is list, should contain estimate and se")
					lambda <- clambda
					lambda$se <- NA
				} else {
					stop("clambda should be logical, numeric, or list")
				}
			}
			chi2.c1df <- chi2.1df/lambda$estimate

			if (is.logical(clambda)) {
				lambda$iz0 <- estlambda(z0*z0,plot=FALSE,prop=propPs)$estimate 
				lambda$iz2 <- estlambda(z2*z2,plot=FALSE,prop=propPs)$estimate
				if (clambda && lambda$iz0<1.0) {warning("z0 lambda < 1, set to 1");lambda$iz0<-1.0}
				if (clambda && lambda$iz2<1.0) {warning("z2 lambda < 1, set to 1");lambda$iz2<-1.0}
				chi2.c2df <- (z0*z0/lambda$iz0 + z2*z2/lambda$iz2 - 2.*z0*z2*rho/(sqrt(lambda$iz0*lambda$iz2)))/(1.- rho*rho)
			} else {
				if (is.list(clambda) && !any(is.na(match(c("estimate","iz0","iz2"),names(clambda))))) {
					chi2.c2df <- (z0*z0/lambda$iz0 + z2*z2/lambda$iz2 - 2.*z0*z2*rho/(sqrt(lambda$iz0*lambda$iz2)))/(1.- rho*rho)
				} else {
					lambda$iz0 <- 1.0
					lambda$iz2 <- 1.0
					chi2.c2df <- chi2.2df
				}
			}
			effB <- chi2[(3*lenn+1):(lenn*4)]
			effB[abs(effB+999.99)<1.e-8] <- NA
			effAB <- chi2[(4*lenn+1):(lenn*5)]
			effAB[abs(effAB+999.99)<1.e-8] <- NA
			effBB <- chi2[(5*lenn+1):(lenn*6)]
			effBB[abs(effBB+999.99)<1.e-8] <- NA
			if (times>1) {
				pr.1df <- rep(0,lenn)
				pr.2df <- rep(0,lenn)
				pr.c1df <- rep(0,lenn)
				pr.c2df <- rep(0,lenn)
			}
		} else {
			th1 <- max(chi2[1:lenn])
			pr.1df <- pr.1df + 1*(chi2.1df < th1)
			pr.2df <- pr.2df + 1*(chi2.2df < max(chi2[(lenn+1):(2*lenn)]))
			pr.c1df <- pr.c1df + 1*(chi2.c1df < th1)
			pr.c2df <- pr.c2df + 1*(chi2.c2df < th1)
			if (!quiet && ((j-1)/bcast == round((j-1)/bcast))) {
#				cat("\b\b\b\b\b\b",round((100*(j-1)/times),digits=2),"%",sep="")
				cat(" ",round((100*(j-1)/times),digits=2),"%",sep="")
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
		out$Pc2df <- pr.c2df/times
		out$Pc2df <- replace(out$Pc2df,(out$Pc2df==0),1/(1+times))
	} else {
		out$P1df <- pchisq(chi2.1df,1,lower=F)
		out$P2df <- pchisq(chi2.2df,actdf,lower=F)
		out$Pc1df <- pchisq(chi2.c1df,1,lower=F)
		out$Pc2df <- pchisq(chi2.c2df,2,lower=F)
	}
	out$lambda <- lambda
	out$effB <- effB
	out$effAB <- effAB
	out$effBB <- effBB
	out$N <- chi2[(6*lenn+1):(lenn*7)]
	if (details) {
		out$snpnames <- data@gtdata@snpnames
		out$idnames <- data@gtdata@idnames
	}
	out$map <- data@gtdata@map
	out$chromosome <- data@gtdata@chromosome
	out$formula <- match.call()
	out$family <- paste("score test for association with trait type",trait.type)
	class(out) <- "scan.gwaa"
	out
}

###### ----------------------- ##################

"qtscore.old" <-
function(formula,data,snpsubset,idsubset,strata,trait.type="gaussian",times=1,quiet=FALSE,bcast=10,clambda=TRUE,propPs=1.0,details=TRUE) {
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
				if (bin == 1) warning("binomial traits is analysed as gaussian")
				bin <- 0
			} else {
				stop("trait type should be either binomial or gaussian")
			}
		}
		if (bin == 1) {
			trait.type="binomial" 
			tmp <- levels(as.factor(y))
			if (tmp[1] != "0" || tmp[2] != "1") stop("binomial outcome should be coded as 0/1")
		} else if (bin == 0) {
			trait.type="gaussian"
		}
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
		chi2 <- .C("qtscore",as.raw(data@gtdata@gtps),as.double(resid),as.integer(bin),as.integer(data@gtdata@nids),as.integer(data@gtdata@nsnps), as.integer(nstra), as.integer(strata), chi2 = double(7*data@gtdata@nsnps), PACKAGE="GenABEL")$chi2
		if (any(data@gtdata@chromosome=="X")) {
		  ogX <- data@gtdata[,data@gtdata@chromosome=="X"]
		  sxstra <- strata; sxstra[ogX@male==1] <- strata[ogX@male==1]+nstra
		  chi2X <- .C("qtscore",as.raw(ogX@gtps),as.double(resid),as.integer(bin),as.integer(ogX@nids),as.integer(ogX@nsnps), as.integer(nstra*2), as.integer(sxstra), chi2 = double(7*ogX@nsnps), PACKAGE="GenABEL")$chi2
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
		out$P1df <- pchisq(chi2.1df,1,lower=F)
		out$P2df <- pchisq(chi2.2df,actdf,lower=F)
		out$Pc1df <- pchisq(chi2.c1df,1,lower=F)
	}
	out$lambda <- lambda
	out$effB <- effB
	out$effAB <- effAB
	out$effBB <- effBB
	out$N <- chi2[(6*lenn+1):(lenn*7)]
	if (details) {
		out$snpnames <- data@gtdata@snpnames
		out$idnames <- data@gtdata@idnames
	}
	out$map <- data@gtdata@map
	out$chromosome <- data@gtdata@chromosome
	out$formula <- match.call()
	out$family <- paste("score test for association with trait type",trait.type)
	class(out) <- "scan.gwaa"
	out
}

"test.type" <-
function(resid,trait.type) {
	if (ismono(resid)) stop("trait is monomorphic")
	if (isbinomial(resid)) {
		if (trait.type == "gaussian") warning("binomial trait is analysed as gaussian")
		tmp <- levels(as.factor(resid))
		if (tmp[1] != "0" || tmp[2] != "1") stop("binomial outcome should be coded as 0/1")
	} else {
		if (trait.type == "binomial") stop("can not analyse trait with > 2 levels as binomial")
	}
}

"isbinomial" <- 
function(y) {
	y <- y[!is.na(y)]
	if (length(unique(y)) > 2) return(FALSE)
	else return(TRUE) 
}

"ismono" <- 
function(y) {
	y <- y[!is.na(y)]
	if (length(unique(y)) <= 1) return(TRUE)
	else return(FALSE) 
}
