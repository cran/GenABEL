"check.marker" <-
function(data, snpsubset, idsubset,
			callrate=0.95,perid.call=0.95, extr.call = 0.1, extr.perid.call = 0.1, 
			het.fdr=0.01, ibs.threshold = 0.95, ibs.mrk = 2000, ibs.exclude="lower",
			maf, p.level=-1, 
			fdrate = 0.2, odds = 1000, hweidsubset, redundant="no", minconcordance = 2.0, 
			qoption="bh95",imphetasmissing=TRUE,XXY.call=0.8) {

	if (class(data) == "gwaa.data") data <- data@gtdata
	if (class(data) != "snp.data") stop("data argument should be of type gwaa.data or snp.data");
	if (!missing(snpsubset)) data <- data[,snpsubset]
	if (!missing(idsubset)) data <- data[idsubset,]
	if (!missing(hweidsubset)) {
		if (is.logical(hweidsubset)) hweidsubset <- data@idnames[which(hweidsubset==TRUE)]
		if (is.numeric(hweidsubset)) hweidsubset <- data@idnames[hweidsubset[!is.na(hweidsubset)]]
	}
	gc(verbose=FALSE)
	
	out <- list()
	snam <- data@snpnames
	sids <- data@idnames
	smap <- data@map
	schr <- data@chromosome

	cat("Excluding people/markers with extremely low call rate...\n")
	cat(data@nsnps,"markers and",data@nids,"people in total\n")
	ts <- perid.summary(data)
	out$idnocall <- rownames(ts[ts[,"CallPP"]<extr.perid.call,])
	out$idok <- data@idnames[!(data@idnames %in% out$idnocall)]
	cat(length(out$idnocall),"people excluded because of call rate <",extr.perid.call,"\n")
	ts <- summary(data)
	if (any(as.character(data@chromosome) == "Y") && any(data@male==1)) {
		tsY <- summary(data[which(data@male==1),(which(as.character(data@chromosome) == "Y"))])
		ts[(which(as.character(data@chromosome) == "Y")),] <- tsY
	}
	out$nocall <- rownames(ts[ts[,"CallRate"]<extr.call,])
	out$snpok <- data@snpnames[!(data@snpnames %in% out$nocall)]
	cat(length(out$nocall),"markers excluded because of call rate <",extr.call,"\n")
	data <- data[out$idok,out$snpok]
	gc(verbose=FALSE)
	cat("Passed:",length(out$snpok),"markers and",length(out$idok),"people\n")
	if (!missing(hweidsubset)) hweidsubset <- hweidsubset[hweidsubset %in% out$idok]
	
	updat <- 0
	if (any(data@chromosome=="X")) {
		cat("\nRunning sex chromosome checks...\n")
		out.nxt <- Xcheck(data[,data@chromosome=="X"],Pgte=0.001,Pssw=0.01,Pmsw=0.01,odds=odds,tabonly=F)
		nxerr <- dim(out.nxt$Xerrtab)[1]
		if (!length(nxerr)) nxerr <- 0
		cat(nxerr,"heterozygous X-linked male genotypes found\n")
		cat(length(out.nxt$Xmrkfail),"X-linked markers are likely to be autosomal (odds >",odds,")\n")
		cat(length(out.nxt$isfemale),"male are likely to be female (odds >",odds,")\n")
		cat(length(out.nxt$ismale),"female are likely to be male (odds >",odds,")\n")
		out <- update.check.marker(out,out.nxt)
		updat <- 1
		out.nxt <- Xcheck(data[out$idok,out$snpok[out$snpok %in% data@snpnames[data@chromosome=="X"]]],Pgte=0.001,Pssw=0.01,Pmsw=0.01,odds=odds,tabonly=T)
		nxerr <- dim(out.nxt$Xerrtab)[1]
		if (!length(nxerr)) nxerr <- 0
		cat("If these people/markers are removed,",nxerr,"heterozygous male genotypes are left\n")
		if (nxerr & imphetasmissing) cat("... these will be considered missing in analysis.\n")
		if (nxerr) cat("... Use Xfix() to fix these problems.\n")
	}

# XXY checks
	if (any(data@chromosome=="Y")) {
		totest <- unique(c(out$isfemale,intersect(data@idnames[data@male==0],out$idok)))
		Ytab <- perid.summary(data[totest,(as.character(data@chromosome) == "Y")])
		Ytab <- Ytab[Ytab$NoMeasured>0,]
		cat("\n",sum(Ytab$NoMeasured),"possibly female Y genotypes identified")
		Ytabexc <- Ytab[Ytab$CallPP>=XXY.call,]
		Ytabnot <- Ytab[Ytab$CallPP<XXY.call,]
		if (dim(Ytabexc)[1]>0) {
			out.nxt$isXXY <- rownames(Ytabexc)
			cat("\n",length(out.nxt$isXXY),"people excluded as possible XXY (Y call >=",XXY.call,")\n")
			out <- update.check.marker(out,out.nxt)
			updat <- 1
			if (dim(Ytabnot)[1]>0) {
				cat("If these people are excluded",sum(Ytabnot$NoMeasured),"female Y-genotypes are left\n")
				if (nxerr & imphetasmissing) cat("... these will be considered missing in analysis.\n")
				if (nxerr) cat("... Use Xfix() to fix these problems.\n")
			}
		} else {
			cat("\nNone of these people excluded based on Y-threshold of",XXY.call,"\n")
		}
	}

	if (updat) { 
		data <- data[out$idok,out$snpok]
		gc(verbose=FALSE)
		cat("Passed:",length(out$snpok),"markers and",length(out$idok),"people\n")
	}
#	out<-list(qc=out,data=data)
#	return(out)

	ts <- summary(data)
	if (any(data@chromosome=="Y")) {
		cat("\nChecking Y-chromsome heterozygous genotypes... ");
		sumYerr <- sum(ts[data@chromosome=="Y","P.12"],na.rm=T)
		cat(sumYerr,"(",100*sumYerr/sum(ts[data@chromosome=="Y","NoMeasured"],na.rm=T),"%) found.\n");
		if (sumYerr & imphetasmissing) cat("... these will be considered missing in analysis.\n")
		if (sumYerr) cat("... Use Xfix() to fix these problems.\n")
	}
	if (any(data@chromosome=="mt")) {
		cat("\nChecking mtDNA heterozygous genotypes... ");
		sumMTerr <- sum(ts[data@chromosome=="mt","P.12"],na.rm=T)
		cat(sumMTerr,"(",100*sumMTerr/sum(ts[data@chromosome=="mt","NoMeasured"],na.rm=T),"%) found.\n");
		if (sumMTerr & imphetasmissing) cat("... these will be considered missing in analysis.\n")
		if (sumMTerr) cat("... Use Xfix() to fix these problems.\n")
	}
	rm(ts)
	gc(verbose=FALSE)

	if (imphetasmissing & (data@nsnps != length(autosomal(data)))) {
		cat("\n")
		data <- Xfix(data)
		gc(verbose=FALSE)
		cat("\n")
	}
	

# first round -- redundancy 
	run <- 1
	passed <- 0
	while (!passed) {
		cat("\nRUN",run,"\n");
		out.nxt <- check.marker.internal(data=data, 
			callrate=callrate,perid.call=perid.call,
			het.fdr=het.fdr, ibs.threshold = ibs.threshold, ibs.mrk = ibs.mrk, ibs.exclude=ibs.exclude,
			maf = maf, p.level= p.level, 
			fdrate = fdrate, hweidsubset = hweidsubset, redundant = "no", minconcordance = 2.0, 
			qoption = qoption, imphetasmissing = imphetasmissing)
		if (length(out$snpok) == length(out.nxt$snpok) && length(out$idok) == length(out.nxt$idok))
			if (all(out$snpok == out.nxt$snpok) && all(out$idok == out.nxt$idok)) 
				passed <- 1
		if (!passed) {
			out <- update.check.marker(out,out.nxt)
			data <- data[out$idok,out$snpok]
			gc(verbose=FALSE)
			if (!missing(hweidsubset)) hweidsubset <- hweidsubset[hweidsubset %in% out$idok]
			run <- run + 1
		}
	}
	if (redundant != "no") {
		cat("\nCHECK REDUNDANCY\n");
		out.nxt <- check.marker.internal(data=data, 
			callrate=callrate,perid.call=perid.call,
			het.fdr=het.fdr, ibs.threshold = ibs.threshold, ibs.mrk = ibs.mrk, ibs.exclude=ibs.exclude,
			maf = maf, p.level= p.level, 
			fdrate = fdrate, hweidsubset = hweidsubset, redundant = redundant, minconcordance = minconcordance, 
			qoption = qoption, imphetasmissing = imphetasmissing)
		out <- update.check.marker(out,out.nxt)
	}
	out$call$name <- snam
	out$call$ids <- sids
	out$call$map <- smap
	out$call$chromosome <- schr
	out$call$call <- match.call()
	class(out) <- "check.marker"
	out
}
