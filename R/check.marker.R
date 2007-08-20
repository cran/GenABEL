"check.marker" <-
function(data, snpsubset, idsubset,
			callrate=0.95,perid.call=0.95, extr.call = 0.1, extr.perid.call = 0.1, 
			het.fdr=0.01, ibs.threshold = 0.95, ibs.mrk = 2000, ibs.exclude="lower",
			maf, p.level=-1, 
			fdrate = 0.2, odds = 1000, hweidsubset, redundant="no", minconcordance = 2.0, 
			qoption="bh95") {

	if (class(data) == "gwaa.data") data <- data@gtdata
	if (class(data) != "snp.data") stop("data argument should be of type gwaa.data or snp.data");
	if (!missing(snpsubset)) data <- data[,snpsubset]
	if (!missing(idsubset)) data <- data[idsubset,]
	if (!missing(hweidsubset)) {
		if (is.logical(hweidsubset)) hweidsubset <- data@idnames[which(hweidsubset==TRUE)]
		if (is.numeric(hweidsubset)) hweidsubset <- data@idnames[hweidsubset[!is.na(hweidsubset)]]
	}
	
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
	out$nocall <- rownames(ts[ts[,"CallRate"]<extr.call,])
	out$snpok <- data@snpnames[!(data@snpnames %in% out$nocall)]
	cat(length(out$nocall),"markers excluded because of call rate <",extr.call,"\n")
	cat("Passed:",length(out$snpok),"markers and",length(out$idok),"people\n")
	data <- data[out$idok,out$snpok]
	if (!missing(hweidsubset)) hweidsubset <- hweidsubset[hweidsubset %in% out$idok]
	
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
		data <- data[out$idok,out$snpok]
		out.nxt <- Xcheck(data[,data@chromosome=="X"],Pgte=0.001,Pssw=0.01,Pmsw=0.01,odds=odds,tabonly=T)
		nxerr <- dim(out.nxt$Xerrtab)[1]
		if (!length(nxerr)) nxerr <- 0
		cat("If these people/markers are removed,",nxerr,"heterozygous male genotypes are left\n")
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
			qoption = qoption)
		if (length(out$snpok) == length(out.nxt$snpok) && length(out$idok) == length(out.nxt$idok))
			if (all(out$snpok == out.nxt$snpok) && all(out$idok == out.nxt$idok)) 
				passed <- 1
		if (!passed) {
			out <- update.check.marker(out,out.nxt)
			data <- data[out$idok,out$snpok]
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
			qoption = qoption)
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
