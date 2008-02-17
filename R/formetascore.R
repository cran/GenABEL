"formetascore" <-
function(formula,data,stat=qtscore,transform=ztransform,build="unknown",verbosity=1, ...) {
	if (is.character(transform)) {
		if (transform=="no") {
			attach(data@phdata,pos=2,warn.conflicts=FALSE)
			if (class(formula) == "formula") {
				mf <- model.frame(formula,data@phdata,na.action=na.omit,drop.unused.levels=TRUE)
				y <- model.response(mf)
				desmat <- model.matrix(formula,mf)
				tmp <- pmatch(names(match.call()),"family")
				tmp <- tmp[!is.na(tmp)]
				if (any(tmp==1)) {
					lmf <- glm.fit(desmat,y,family=family)
				} else {
					lmf <- glm.fit(desmat,y,family=gaussian())
				}
				mids <- rownames(data) %in% rownames(mf)
				resid <- lmf$resid
			} else if (class(formula) == "numeric" || class(formula) == "integer" || class(formula) == "double") {
				y <- formula
				mids <- (!is.na(y))
				y <- y[mids]
				resid <- y
				if (length(unique(resid))==1) stop("trait is monomorphic")
				if (length(unique(resid))==2) stop("trait is binary")
			} else {
				stop("formula argument must be a formula or one of (numeric, integer, double)")
			}
			detach(data@phdata)
			trvar <- resid
		} else {
			stop("transform argument not recognised")
		}
	} else {
		trvar <- transform(formula=formula,data=data)
	}
	if (verbosity<0) stop("verbosity parameter must be positive integer")
	data <- data[!is.na(trvar),]
	trvar <- trvar[!is.na(trvar)]
	res <- stat(trvar,data,...)
	sum <- summary(data@gtdata)
	callr <- sum$Call
	Phwe <- sum$Pexact
	efff <- sum$Q.2
	reff <- 1. - sum$Q.2
	serr <- abs(res$effB)/sqrt(res$chi2.1df)
	coding <- as.character(data@gtdata@coding)
	if (verbosity==0)
 	out <- data.frame(name=res$snpnames,chromosome=res$chromosome,
		position=res$map,strand=as.character(data@gtdata@strand),
		allele1=alleleID.reference()[coding],
		allele2=alleleID.effective()[coding],
		effallele=alleleID.effective()[coding],
		beta=res$effB,sebeta=serr,p=res$P1df,
		stringsAsFactors=F)
	else if (verbosity==1)
 	out <- data.frame(name=res$snpnames,chromosome=res$chromosome,
		position=res$map,strand=as.character(data@gtdata@strand),
		allele1=alleleID.reference()[coding],
		allele2=alleleID.effective()[coding],
		effallele=alleleID.effective()[coding],
		effallelefreq=efff,
		n=res$N,beta=res$effB,sebeta=serr,p=res$P1df,
		pgc=res$Pc1df,
		pexhwe=Phwe,call=callr,stringsAsFactors=F)
	else 
 	out <- data.frame(name=res$snpnames,chromosome=res$chromosome,
		position=res$map,strand=as.character(data@gtdata@strand),
		allele1=alleleID.reference()[coding],
		allele2=alleleID.effective()[coding],
		build=rep(build,length(coding)),
		effallele=alleleID.effective()[coding],
		effallelefreq=efff,
		n=res$N,beta=res$effB,sebeta=serr,p=res$P1df,
		pgc=res$Pc1df,lambda=rep(res$lam$est,length(coding)),
		pexhwe=Phwe,call=callr,fmax=sum$Fmax,plrthwe=sum$Plrt,stringsAsFactors=F)

	rownames(out) <- res$snpnames
	out
}
