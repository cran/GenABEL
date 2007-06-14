"polygenic" <-
function(formula,kinship.matrix,data,fixh2,starth2=0.3,trait.type="gaussian",opt.method="nlm",scaleh2=1000,quiet=FALSE,...) {
	if (!missing(data)) if (class(data) == "gwaa.data") data <- data@phdata
	if (!missing(data)) if (class(data) != "data.frame") stop("data should be of gwaa.data or data.frame class")
	if (starth2<0 || starth2>1) stop("Starting value of h2 (starth2) must be between 0 and 1")
	ttargs <- c("gaussian","binomial")
	if (!(match(trait.type,ttargs,nomatch=0)>0)) {
		out <- paste("trait.type argument should be one of",ttargs,"\n")
		stop(out)
	}
	if (trait.type=="gaussian") fam <- gaussian()
	if (trait.type=="binomial") fam <- binomial()
	optargs <- c("nlm","optim")
	if (!(match(opt.method,optargs,nomatch=0)>0)) {
		out <- paste("opt.method argument should be one of",optargs,"\n")
		stop(out)
	}
	
	if (!missing(data)) attach(data,pos=2,warn.conflicts=FALSE)
	if (class(formula) == "formula") {
		clafo <- "formula"
		mf <- model.frame(formula,data,na.action=na.omit,drop.unused.levels=TRUE)
		y <- model.response(mf)
		sdy <- sd(y)
		meany <- mean(y)
		if (trait.type=="gaussian") y <- (y-meany)/sdy
		desmat <- model.matrix(formula,mf)
		lmf <- glm.fit(desmat,y,family=fam)
		iniest <- lmf$coeff
		mids <- rownames(data) %in% rownames(mf)
		if (trait.type=="binomial") {
			tvar <- var(lmf$resid)
		} else {
			tvar <- var(lmf$resid)
		}
	} else {
		clafo <- "NOT"
		y <- formula
		mids <- (!is.na(y))
		y <- y[mids]
		sdy <- sd(y)
		meany <- mean(y)
		if (trait.type=="gaussian") y <- (y-meany)/sdy
		desmat <- matrix(1,nrow=length(y))
		if (trait.type=="binomial") {
			tvar <- var(y)
			tmp <- mean(y)
			iniest <- c(log(tmp/(1.-tmp)))
		} else {
			tvar <- var(y)
			iniest <- c(mean(y))
		}
	}
	if (!missing(data)) detach(data)
	relmat <- kinship.matrix[mids,mids]*2.0
	tmp <- t(relmat)
	relmat[upper.tri(relmat)] <- tmp[upper.tri(tmp)]
	rm(tmp);gc()
	eigres <- eigen(ginv(relmat),symm=TRUE)
	print(iniest);flush.console()
	if (!missing(fixh2)) {
		startlik <- polylik(c(iniest,tvar),y=y,desmat=desmat,relmat=relmat,ervec=eigres$vec,fixh2=(fixh2/scaleh2),trait.type=trait.type,scaleh2=scaleh2)
		if (opt.method=="nlm") {
		prnlev <- 0; if (!quiet) prnlev <- 2;
		h2an <- nlm(polylik,p=c(iniest,tvar),y=y,desmat=desmat,relmat=relmat,ervec=eigres$vec,fixh2=(fixh2/scaleh2),trait.type=trait.type,startlik=startlik,scaleh2=scaleh2,print.level=prnlev,...)
		} else {
		lower <- c(rep(-scaleh2,length(iniest)),1.e-16)
		upper <- c(rep(scaleh2,length(iniest)),(2*scaleh2))
		par <- c(iniest,starth2,tvar)
		cntrl <- list(); if (!quiet) cntrl <- list(trace=6,REPORT=1)
		h2an <- optim(fn=polylik,par=c(iniest,tvar),method="L-BFGS-B",lower=c(rep(-Inf,length(iniest)),1.e-16),upper=c(rep(Inf,length(iniest)),Inf),y=y,desmat=desmat,relmat=relmat,ervec=eigres$vec,fixh2=(fixh2),trait.type=trait.type,control=cntrl,scaleh2=1,...)
		}
	} else {
		startlik<-polylik(c(iniest,(starth2/scaleh2),tvar),y=y,desmat=desmat,relmat=relmat,ervec=eigres$vec,trait.type=trait.type,scaleh2=scaleh2)
		if (opt.method=="nlm") {
		prnlev <- 0; if (!quiet) prnlev <- 2;
		h2an <- nlm(polylik,p=c(iniest,(starth2/scaleh2),tvar),y=y,desmat=desmat,relmat=relmat,ervec=eigres$vec,trait.type=trait.type,startlik=startlik,scaleh2=scaleh2,print.level=prnlev,...)
		} else {
		lower <- c(rep(-scaleh2,length(iniest)),1.e-16,1.e-16)
		upper <- c(rep(scaleh2,length(iniest)),1.-1.e-16,(2*scaleh2))
		par <- c(iniest,starth2,tvar)
		cntrl <- list(); if (!quiet) cntrl <- list(trace=6,REPORT=1)
		h2an <- optim(fn=polylik,par=par,method="L-BFGS-B",lower=lower,upper=upper,y=y,desmat=desmat,relmat=relmat,ervec=eigres$vec,trait.type=trait.type,control=cntrl,scaleh2=1,...)
		}
	}
	if (opt.method=="optim") {
		scaleh2 <- 1
		h2an$estimate <- h2an$par; h2an$par
		h2an$minimum <- h2an$value; h2an$value
	}
	out <- list();
	npar <- length(h2an$estimate)
	if (trait.type=="gaussian") {
		if (!missing(fixh2)) {
			h2an$estimate[1:(npar-1)] <- h2an$estimate[1:(npar-1)]*sdy
			h2an$estimate[npar] <- h2an$estimate[npar]*sdy*sdy
			h2an$estimate[1] <- h2an$estimate[1]+meany
		} else {
			h2an$estimate[1:(npar-2)] <- h2an$estimate[1:(npar-2)]*sdy
			h2an$estimate[(npar-1)] <- h2an$estimate[(npar-1)]*scaleh2
			h2an$estimate[npar] <- h2an$estimate[npar]*sdy*sdy
			h2an$estimate[1] <- h2an$estimate[1]+meany
		}
	}
	out$h2an <- h2an
	if (trait.type=="gaussian") {scay <- y*sdy+meany} else {scay<-y}
	if (clafo == "formula") {
		if (!missing(fixh2)) {fxeff <- h2an$est[1:(npar-1)]} else {fxeff <- h2an$est[1:(npar-2)]}
		if (trait.type=="gaussian") {eY <- desmat %*% fxeff} else {ee <- exp(desmat %*% fxeff); eY <- ee/(1.+ee);}
		out$residualY <- scay - eY
	} else {
		if (trait.type=="gaussian") {eY <- h2an$estimate[1]} else {ee <- exp(h2an$estimate[1]); eY <- ee/(1.+ee);}
		out$residualY <- scay - eY
	}
	if (!missing(fixh2)) {
		h2 <- fixh2
	} else {
		h2 <- h2an$estimate[npar-1]
	}
	out$esth2 <- h2
	tvar <- h2an$estimate[npar]
	ervec <- eigres$vec
	sigma <- h2*tvar*relmat + (1-h2)*tvar*diag(x=1,ncol=length(y),nrow=length(y))
#	es <- 1./diag(t(ervec) %*% (sigma) %*% ervec)
#	ginvsig <- ervec %*% diag(es,ncol=length(y)) %*% t(ervec)
	out$InvSigma <- ginv(sigma) #ginvsig
	out$pgresidualY <- as.vector((1.-h2) * tvar * (out$InvSigma %*% out$residualY))
	out$call <- match.call()
	out$measuredIDs <- mids
	out
}

"polylik" <- function(coeffs,y,desmat,relmat,ervec,fixh2,trait.type,startlik=.Machine$double.xmax,scaleh2=1) {
	if (!missing(fixh2)) {
		h2 <- fixh2
		fixeff <- coeffs[1:(length(coeffs)-1)]
	} else {
		h2 <- coeffs[(length(coeffs)-1)] # coeffs[(length(coeffs))] 
		fixeff <- coeffs[1:(length(coeffs)-2)]
	}
	h2 <- h2*scaleh2
	tvar <- coeffs[(length(coeffs))] # var(y) 
	if (h2<0 || h2>1 || tvar<0) return(startlik)
	nids <- length(y)
	sigma <- h2*tvar*relmat + (1-h2)*tvar*diag(x=1,ncol=nids,nrow=nids)
	if (trait.type=="binomial") {
		expb <- exp(desmat %*% fixeff)
		qt <- y - expb/(1.+expb)
	} else {
		qt <- y - desmat %*% fixeff
	}
	es <- 1./diag(t(ervec) %*% (sigma) %*% ervec)
	ginvsig <- ervec %*% diag(es,ncol=length(qt)) %*% t(ervec)
	a <- determinant(sigma,logarithm=T)
	a <- a$modulus * a$sign
	b <- (t(qt) %*% ginvsig %*% qt)
	loglik <- a+b
	loglik
}

