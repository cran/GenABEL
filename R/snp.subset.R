snp.subset <- function(data,snpsubset) {
	if (missing(snpsubset)) {
		warning("snpsubset missing")
		return(data)
	}
	if (class(data) == "scan.gwaa") {
		out <- snp.subset.scan.gwaa(data,snpsubset)
	} else if (class(data) == "check.marker") {
		out <- snp.subset.check.marker(data,snpsubset)
	} else {stop("data should be of type check.marker-class or scan.gwaa-class")}
	out
}

"snp.subset.scan.gwaa" <- 
function(data,snpsubset) {
	norig <- length(data$snpnames)
	tokeep <- rep(FALSE,norig)
	if (is.numeric(snpsubset) || is.logical(snpsubset)) {
		tokeep[snpsubset] <- TRUE
	} else if (is.character(snpsubset)) {
		tokeep <- data$snpnames %in% snpsubset
	} else {stop("snpsubset should have type numeric, logical or character")}
	out <- list(P1df=data$P1df[tokeep],P2df=data$P2df[tokeep],Pc1df=data$Pc1df[tokeep],chi2.1df=data$chi2.1df[tokeep],chi2.2df=data$chi2.2df[tokeep],N=data$N[tokeep], snpnames=data$snpnames[tokeep],formula=match.call(),family=data$family,map=data$map[tokeep],idnames=data$idnames,chromosome=data$chromosome[tokeep],lambda=data$lambda,effB=data$effB[tokeep],effAB=data$effAB[tokeep],effBB=data$effBB[tokeep])
	class(out) <- "scan.gwaa"
	out 
}

"snp.subset.check.marker" <- 
function(mcobj,subs) {
	o <- list()
# processing mcobj$call
	cidx <- match(subs,mcobj$call$name)
	o$call$name <- mcobj$call$name[cidx]
	o$call$map <- mcobj$call$map[cidx]
	o$call$chromosome <- mcobj$call$chromosome[cidx]
	o$call$call <- match.call()
# processing mcobj$details.redundancy
	for (i in names(mcobj$details.redundancy)) {
		if (any(subs==i))
		o$details.redundancy[i] <- mcobj$details.redundancy[i]
	}

# processing simple lists
	o$redundant <- subs[!is.na(match(subs,mcobj$redundant))]
	o$nofreq <- subs[!is.na(match(subs,mcobj$nofreq))]
	o$nocall <- subs[!is.na(match(subs,mcobj$nocall))]
	o$ok <- subs[!is.na(match(subs,mcobj$ok))]

# out hwe & chi2
	cidx <- !is.na(match(mcobj$nohwe,subs))
	o$nohwe <- mcobj$nohwe[cidx]
	o$chi2.hwe <- mcobj$chi2.hwe[cidx]

# output
	class(o) <- "check.marker"
	o
}
