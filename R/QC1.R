doQC <- function(formula,data,plot=T,update=T) {
	if (missing(data)) stop("data argument missing")
	if (is.character(formula)) {
		attach(data)
		var <- get(formula);
		detach(data)
	} else {
		lmr <- lm(formula,data=data);
		var <- lmr$res + lmr$coeff["(Intercept)"];
	}
	if (any(var<=1.e-16,na.rm=T)) {
		warning(paste("Results of",formula,"contain zeros and/or negative values;",floor(-min(var)+1),"added to var"));
		var <- var +floor(- min(var) + 1)
	}
	lnvar <- log(var)
	sqvar <- sqrt(var)
	Zvar <- rank(var)-.5
	Zvar <- replace(Zvar,is.na(var),NA)
	mP <- .5/max(Zvar,na.rm=T)
	Zvar <- Zvar/(max(Zvar,na.rm=T)+.5)
	Zvar <- qnorm(Zvar)
	ZZvar <- Zvar
	pred <- ppoints(Zvar)
	ZZvar[!is.na(Zvar)] <- lm(Zvar~pred)$res
	if (plot==T) {
		par(mfcol=c(2,5))
		hist(var,main=as.character(formula));
		qqnorm(var); qqline(var)
		hist(lnvar,main=paste("LOG",as.character(formula)));
		qqnorm(lnvar); qqline(lnvar)
		hist(sqvar,main=paste("SQRT",as.character(formula)));
		qqnorm(sqvar); qqline(sqvar)
		hist(Zvar,main=paste("Z-trans",as.character(formula)));
		qqnorm(Zvar); qqline(Zvar)
		hist(ZZvar,main=paste("ZZ-trans",as.character(formula)));
		qqnorm(ZZvar); qqline(ZZvar)
	}
	x<-(shapiro.test(var))
	cat(paste("Raw",as.character(formula),": W =",round(x$stat,dig=5),"; P =",x$p),"\n")
	x<-(shapiro.test(lnvar))
	cat(paste("LOG",as.character(formula),": W =",round(x$stat,dig=5),"; P =",x$p),"\n")
	x<-(shapiro.test(sqvar))
	cat(paste("SQRT",as.character(formula),": W =",round(x$stat,dig=5),"; P =",x$p),"\n")
	x<-(shapiro.test(Zvar))
	cat(paste("Z-TR",as.character(formula),": W =",round(x$stat,dig=5),"; P =",x$p),"\n")
	x<-(shapiro.test(ZZvar))
	cat(paste("ZZ-TR",as.character(formula),": W =",round(x$stat,dig=5),"; P =",x$p),"\n")
	out <-list()
	out$var <- var;
	out$sqvar <- sqvar;
	out$lnvar <- lnvar;
	out$Zvar <- Zvar;
	if (is.character(formula) && !missing(data) && update==T) {
		data$tmpsq <- sqvar;
		data$tmpln <- lnvar;
		data$tmpZ <- Zvar;
		data$tmpZZ <- ZZvar;
		x <- names(data)
		x[(length(x)-3)] <- paste("sq",formula,sep="")
		x[(length(x)-2)] <- paste("ln",formula,sep="")
		x[(length(x)-1)] <- paste("ztold",formula,sep="")
		x[length(x)] <- paste("zt",formula,sep="")
		names(data) <- x
	}
	data
}


