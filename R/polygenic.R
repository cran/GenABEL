#' Estimation of polygenic model
#' 
#'  This function maximises the likelihood of the data under polygenic 
#' 	model with covariates an reports twice negative maximum likelihood estimates 
#' 	and the inverse of variance-covariance matrix at the point of ML. 
#' 
#' 	One of the major use of this function is to estimate residuals of the 
#' 	trait and the inverse of the variance-covariance matrix for 
#' 	further use in analysis with \code{\link{mmscore}} and 
#' 	\code{\link{grammar}}.
#' 
#' 	Also, it can be used for a variant of GRAMMAR analysis, which 
#' 	allows for permutations for GW significance by use of 
#' 	environmental residuals as an analysis trait with \code{\link{qtscore}}.
#' 
#' 	"Environmental residuals" (not to be mistaken with just "residuals") are 
#' 	the residual where both the effect of covariates AND the estimated 
#' 	polygenic effect (breeding values) are factored out. This thus 
#' 	provides an estimate of the trait value contributed by environment
#' 	(or, turning this other way around, the part of trait not explained 
#' 	by covariates and by the polygene). Polygenic residuals are estimated 
#' 	as
#' 
#' 	\deqn{
#' 	\sigma^2  V^{-1} (Y - (\hat{\mu} + \hat{\beta} C_1 + ...))
#' 	}
#' 
#' 	where \eqn{sigma^2} is the residual variance, \eqn{V^{-1}} is the 
#' 	InvSigma (inverse of the var-cov matrix at the maximum of 
#' 	polygenic model) and 
#' 	\eqn{(Y - (\hat{\mu} + \hat{\beta} C_1 + ...))} is the trait 
#' 	values adjusted for covariates (also at at the maximum of 
#' 	polygenic model likelihood). 
#' 
#' 	It can also be used for heritability analysis.
#' 	If you want to test significance of heritability, 
#' 	estimate the model and write down 
#' 	the function minimum reported at "h2an" element of the output 
#' 	(this is twice negative MaxLikleihood). Then do next round of 
#' 	estimation, but set fixh2=0. The difference between you function minima 
#' 	gives a test distribued as chi-squared with 1 d.f.
#' 
#' 	The way to compute the likleihood is partly based on 
#' 	the paper of Thompson (see refs), namely instead of 
#' 	taking inverse of var-cov matrix every time, 
#' 	eigenvectors of the inverse of G (taken only once) 
#' 	are used.
#' 
#' 
#'  @param formula Formula describing fixed effects to be used in analysis, e.g. 
#' 	y ~ a + b means that outcome (y) depends on two covariates, a and b. 
#' 	If no covariates used in analysis, skip the right-hand side of the 
#' 	equation.
#'  @param kinship.matrix Kinship matrix, as provided by e.g. ibs(,weight="freq"), 
#' 	or estimated outside of GenABEL from pedigree data.
#'  @param data An (optional) object of \code{\link{gwaa.data-class}} or a data frame with 
#' 	outcome and covariates
#'  @param fixh2 Optional value of heritability to be used, instead of maximisation. 
#' 	The uses of this option are two-fold: (a) testing significance of 
#' 	heritability and (b) using a priori known heritability to derive the 
#' 	rest of MLEs and var.-cov. matrix.
#'  @param starth2 Starting value for h2 estimate
#'  @param trait.type "gaussian" or "binomial"
#'  @param opt.method "nlm" or "optim". These two use dirrerent optimisation functions. 
#'  We suggest using the default \code{\link{nlm}}, though 
#' 	\code{\link{optim}} may give better results in some situations
#'  @param scaleh2 Only relevant when "nlm" optimisation function is used. 
#' 	"scaleh2" is the heritability 
#' 	scaling parameter, regulating how "big" are parameter changes in h2 with the 
#' 	respect to changes in other parameters. As other parameters are estimated 
#' 	from previous regression, these are expected to change little from the 
#' 	initial estimate. The default value of 1000 proved to work rather well under a 
#' 	range of conditions.
#'  @param quiet If FALSE (default), details of optimisation process are reported.
#'  @param steptol steptal parameter of "nlm"
#'  @param gradtol gradtol parameter of "nlm" 
#'  @param optimbou fixed effects boundary scale parameter for 'optim'
#'  @param fglschecks additional check for convergance on/off (convergence 
#'  between estimates obtained and that from FGLS)
#'  @param maxnfgls number of fgls checks to perform
#'  @param maxdiffgls max difference allowed in fgls checks 
#'  @param ... Optional arguments to be passed to \code{\link{nlm}} or (\code{\link{optim}}) 
#' 	minimisation function
#' 
#'  @return 
#'   A list with values 
#'   \item{h2an}{A list supplied by the \code{\link{nlm}} minimisation routine. 
#' 	Of particular interest are elements "estimate" containing parameter 
#' 	maximal likelihood estimates (MLEs) (order: mean, betas for covariates, 
#' 	heritability, (polygenic + residual variance)). The value of 
#' 	twice negative maximum log-likelihood
#' 	is returned as h2an\$minimum.}
#'   \item{residualY}{Residuals from analysis, based on covariate effects only; 
#' 	NOTE: these are NOT grammar "environmental residuals"!}
#'   \item{esth2}{Estimate (or fixed value) of heritability}
#'   \item{pgresidualY}{Environmental residuals from analysis, based on covariate effects 
#' 	and predicted breeding value.
#' 	}
#'   \item{InvSigma}{Inverse of the variance-covariance matrix, computed at the 
#' 	MLEs -- these are used in \code{\link{mmscore}} and \code{\link{grammar}}
#' 	functions.}
#'   \item{call}{The details of call}
#'   \item{measuredIDs}{Logical values for IDs who were used in analysis 
#' 	(traits and all covariates measured) == TRUE}
#'   \item{convFGLS}{was convergence achieved according to FGLS criterionE}
#' 
#' 
#' @references 
#' Thompson EA, Shaw RG (1990) Pedigree analysis for quantitative 
#' traits: variance components without matrix inversion. Biometrics 
#' 46, 399-413.
#' 
#' Aulchenko YS, de Koning DJ, Haley C. Genomewide rapid association using mixed model 
#' and regression: a fast and simple method for genome-wide pedigree-based quantitative 
#' trait loci association analysis. Genetics. 2007 177(1):577-85.
#' 
#' Amin N, van Duijn CM, Aulchenko YS. A genomic background based method for 
#' association analysis in related individuals. PLoS ONE. 2007 Dec 5;2(12):e1274.
#'  
#' @author Yurii Aulchenko
#' 
#' @note 
#' 	Presence of twins may complicate your analysis. Check kinship matrix for 
#' 	singularities, or rather use \code{\link{check.marker}} for identification 
#' 	of twin samples. Take special care in interpretation.
#' 
#' 	If a trait (no covarites) is used, make sure that order of IDs in 
#' 	kinship.matrix is exactly the same as in the outcome
#' 
#' @seealso 
#' \code{\link{mmscore}},
#' \code{\link{grammar}}
#' 
#' @examples 
#' # note that procedure runs on CLEAN data
#' data(ge03d2ex.clean)
#' gkin <- ibs(ge03d2ex.clean,w="freq")
#' h2ht <- polygenic(height ~ sex + age,kin=gkin,ge03d2ex.clean)
#' # estimate of heritability
#' h2ht$esth2
#' # other parameters
#' h2ht$h2an
#' # the minimum twice negative log-likelihood
#' h2ht$h2an$minimum
#' # twice maximum log-likelihood
#' -h2ht$h2an$minimum
#' 
#' #for binary trait (experimental)
#' h2dm <- polygenic(dm2 ~ sex + age,kin=gkin,ge03d2ex.clean,trait="binomial")
#' # estimated parameters
#' h2dm$h2an
#' 
#' @keywords 
#' htest
#' 
"polygenic" <-
		function(formula,kinship.matrix,data,fixh2,starth2=0.3,trait.type="gaussian",
				opt.method="nlm",scaleh2=1,quiet=FALSE,
				steptol=1e-8, gradtol = 1e-8, optimbou = 8, 
				fglschecks=TRUE,maxnfgls=8,maxdiffgls=1e-4, ...) {
	if (!missing(data)) if (is(data,"gwaa.data")) 
		{
			checkphengen(data)
			data <- phdata(data)
		}
	if (!missing(data)) if (!is(data,"data.frame")) stop("data should be of gwaa.data or data.frame class")
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
	if (is(formula,"formula")) {
		print("b")
		clafo <- "formula"
		mf <- model.frame(formula,data,na.action=na.omit,drop.unused.levels=TRUE)
		y <- model.response(mf)
		sdy <- sd(y)
		meany <- mean(y)
		if (trait.type=="gaussian") y <- (y-meany)/sdy
		desmat <- model.matrix(formula,mf)
		rglm <- glm(y~0+desmat,family=fam)
		sglm <- summary(rglm)
		iniest <- sglm$coef[,1]
		inierr <- sglm$coef[,2]
		phids <- rownames(data)[rownames(data) %in% rownames(mf)]
		#print("bb")
		#print(phids)
		relmat <- kinship.matrix[phids,phids]*2.0
		#print("bbb")
		mids <- (rownames(data) %in% rownames(mf))
		if (trait.type=="binomial") {
			tvar <- var(rglm$resid)
		} else {
			tvar <- var(rglm$resid)
		}
	} else {
		clafo <- "NOT"
		y <- formula
		if (length(y) != dim(kinship.matrix)[1]) stop("dimension of outcome and kinship.matrix do not match")
		mids <- (!is.na(y))
		phids <- data$id[mids]
		y <- y[mids]
		relmat <- kinship.matrix[mids,mids]*2.0
		sdy <- sd(y)
		meany <- mean(y)
		if (trait.type=="gaussian") y <- (y-meany)/sdy
		desmat <- matrix(1,nrow=length(y))
		if (trait.type=="binomial") {
			tvar <- var(y)
			tmp <- mean(y)
			iniest <- c(log(tmp/(1.-tmp)))
			inierr <- sqrt(tmp*(1-tmp)/length(y))
		} else {
			tvar <- var(y)
			iniest <- c(mean(y))
			inierr <- sqrt(var(y)/length(y))
		}
	}

	if (!missing(data)) detach(data)
	tmp <- t(relmat)
	relmat[upper.tri(relmat)] <- tmp[upper.tri(tmp)]
	rm(tmp);gc()
	eigres <- eigen(ginv(relmat),symm=TRUE)
	if (!quiet) {
		cat("LM estimates of fixed parameters:\n")
		print(iniest);
		flush.console()
	}
	iniest <- iniest# + 0.001*iniest
	
	
	convFGLS <- NULL;
	
	if (!missing(fixh2)) {
		startlik <- polylik(c(iniest,tvar),y=y,desmat=desmat,relmat=relmat,
				ervec=eigres$vec,fixh2=(fixh2/scaleh2),trait.type=trait.type,
				scaleh2=scaleh2)
		if (opt.method=="nlm") {
			prnlev <- 0; if (!quiet) prnlev <- 2;
			h2an <- nlm(polylik,p=c(iniest,tvar),y=y,desmat=desmat,relmat=relmat,ervec=eigres$vec,
					fixh2=(fixh2/scaleh2),trait.type=trait.type,startlik=startlik,scaleh2=scaleh2,
					print.level=prnlev,steptol=steptol,gradtol=gradtol,...)
#		h2an <- nlm(polylik,p=c(iniest,tvar),y=y,desmat=desmat,relmat=relmat,ervec=eigres$vec,fixh2=(fixh2/scaleh2),trait.type=trait.type,startlik=startlik,scaleh2=scaleh2,...)
		} else {
			lower <- c(iniest-inierr*optimbou,1.e-4)
			upper <- c(iniest+inierr*optimbou,1)
			cntrl <- list(); if (!quiet) cntrl <- list(trace=6,REPORT=1)
			h2an <- optim(fn=polylik,par=c(iniest,tvar),method="L-BFGS-B",lower=lower,upper=upper,y=y,desmat=desmat,relmat=relmat,ervec=eigres$vec,fixh2=(fixh2),trait.type=trait.type,control=cntrl,scaleh2=1,...)
		}
	} else {
		
		
		nfgls <- 0
		diffgls <- maxdiffgls + 1
		parsave <- NULL	
		while (nfgls < maxnfgls & diffgls > maxdiffgls)
		{
			
			
			if (is.null(parsave)) {
				if (opt.method=="nlm") parsave <- c(iniest,starth2/scaleh2,tvar)
				else parsave <- c(iniest,starth2,tvar)
			}
			startlik<-polylik(parsave,y=y,desmat=desmat,relmat=relmat,ervec=eigres$vec,
					trait.type=trait.type,scaleh2=scaleh2)
			if (opt.method=="nlm") {
				#print(parsave)
				prnlev <- 0; if (!quiet) prnlev <- 2;
				h2an <- nlm(polylik,p=parsave,y=y,desmat=desmat,relmat=relmat,ervec=eigres$vec,
						trait.type=trait.type,startlik=startlik,scaleh2=scaleh2,
						print.level=prnlev,steptol=steptol,gradtol=gradtol,...)
#		h2an <- nlm(polylik,p=c(iniest,(starth2/scaleh2),tvar),y=y,desmat=desmat,relmat=relmat,ervec=eigres$vec,trait.type=trait.type,startlik=startlik,scaleh2=scaleh2,...)
				parsave <- h2an$estimate
			} else {
				#print(parsave)
				lower <- c(iniest-inierr*optimbou,1.e-4,1.e-4)
				upper <- c(iniest+inierr*optimbou,0.98,1)
				#print("LOWER/UPPER")
				#print(lower)
				#print(upper)
				cntrl <- list(); if (!quiet) cntrl <- list(trace=6,REPORT=1)
				h2an <- optim(fn=polylik,par=parsave,method="L-BFGS-B",lower=lower,upper=upper,y=y,desmat=desmat,relmat=relmat,ervec=eigres$vec,trait.type=trait.type,control=cntrl,scaleh2=1,...)
				parsave <- h2an$par
			}
			#print(parsave)
			
			
			if (fglschecks && missing(fixh2)) {
				npar <- length(parsave)
				h2 <- parsave[npar-1]*scaleh2
				iSigma <- ginv(h2*relmat + (1-h2)*diag(x=1,ncol=length(y),nrow=length(y)))
				LHest <- parsave[1:(npar-2)]
				betaFGLS <- as.vector(ginv(t(desmat) %*% iSigma %*% desmat) %*% 
								(t(desmat) %*% iSigma %*% y))
				#print(c("betaFGLS = ",betaFGLS))
				#cnd <- (abs(betaFGLS)>maxdiffgls)
				#cnd[1] <- TRUE
				#LHest <- LHest[cnd]
				#betaFGLS <- betaFGLS[cnd]
				difFGLS <- abs(betaFGLS-LHest)
				if (!quiet) {
					cat("difFGLS:\n")
					print(difFGLS)
				}
				diffgls <- max(difFGLS)
				steptol <- steptol/10;
				gradtol <- gradtol/10;
				if ((h2<1e-4 || h2>(1-1e-4)) && nfgls > floor(maxnfgls/2)) {
					parsave[npar-1] <- runif(1,min=0.05,max=0.95)/scaleh2
				}
			} else {
				nfgls <- maxnfgls+1	
			}
			nfgls <- nfgls + 1
			#print(c("NFGLS=",nfgls))
			
		}
		
		if (diffgls<maxdiffgls) {
			convFGLS <- TRUE; 
			if (!quiet) {
				cat("\n******************************************\n");
				cat("*** GOOD convergence indicated by FGLS ***\n");
				cat("******************************************\n");
			}
		} else {
			convFGLS <- FALSE; 
			if (!quiet) {
				cat("\n***********************************************\n");
				cat("*** !!!BAD!!! convergence indicated by FGLS ***\n");
				cat("***********************************************\n");
			}
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
	} else {
		if (!missing(fixh2)) {
		} else {
			h2an$estimate[(npar-1)] <- h2an$estimate[(npar-1)]*scaleh2
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
	rownames(out$InvSigma) <- phids
	colnames(out$InvSigma) <- phids
	pgres <- as.vector((1.-h2) * tvar * (out$InvSigma %*% out$residualY))
	out$measuredIDs <- mids
	names(out$measuredIDs) <- phids
	out$pgresidualY <- rep(NA,length(mids))
	out$pgresidualY[mids] <- pgres
	names(out$pgresidualY) <- phids
	resY <- out$residualY
	out$residualY <- rep(NA,length(mids))
	out$residualY[mids] <- resY
	names(out$residualY) <- phids
	out$call <- match.call()
	out$convFGLS <- convFGLS
	
	a <- determinant(sigma,logarithm=T)
	a <- a$modulus * a$sign
	b <- (t(resY) %*% out$InvSigma %*% resY)
	loglik <- a+b
	out$h2an$minimum <- as.vector(loglik)
	
	
	class(out) <- "polygenic"
	out
}

