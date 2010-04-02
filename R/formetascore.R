"formetascore" <-
function(formula,data,stat=qtscore,transform="no",build="unknown",verbosity=1, ...) {
	if (!is(data,"gwaa.data")) stop("data argument must have gwaa.data-class")
	checkphengen(data)
	if (!missing(data)) attach(data@phdata,pos=2,warn.conflicts=FALSE)
	if (is(formula,"polygenic")) {
		pm <- pmatch("stat",names(match.call()))
		pm <- (pm[!is.na(pm)])[1]
		a <- match.call()[[pm]]
		if (a != "mmscore" && a != "mmscore()") stop("stat should be mmscore when polygenic object is analysed")
		if (!is(transform,"character")) stop("transform should be \"no\" when polygenic object is analysed")
		if (transform != "no") stop("transform should be \"no\" when polygenic object is analysed")
	}
	if (is.character(transform)) {
		if (transform!="no") stop("transform argument not recognised")
	} else {
		formula <- transform(formula=formula,data=data)
	}
	if (is(formula,"formula")) {
		mf <- model.frame(formula,data@phdata,na.action=na.omit,drop.unused.levels=TRUE)
		mids <- rownames(data@phdata) %in% rownames(mf)
	} else if (is(formula,"polygenic")) {
		mids <- which(!is.na(formula$residualY))
	} else {
		mids <- which(!is.na(formula))
	}
	if (!missing(data)) detach(data@phdata)
	if (verbosity<0) stop("verbosity parameter must be positive integer")
	res <- stat(formula,data,...)
	sum <- summary(data@gtdata[mids,])
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
