"check.marker" <-
function(data, snpsubset, idsubset,
			callrate=0.95,maf=0.01, p.level=-1, fdrate = 0.2, hweidsubset, 
			redundant="bychrom", minconcordance = 2.0, 
			qoption="bh95", bcast=10000) {

# qoption = "bh95" (Benjamini & Hochberg 1995) or "storey" (Storey 2003, requires qvalue library)
	if (class(data) == "gwaa.data") {
		if (missing(snpsubset)) snpsubset=data@gtdata@snpnames
		if (missing(idsubset))  idsubset=data@gtdata@idnames
		n1data <- data@gtdata[idsubset,snpsubset]
	} else if (class(data) == "snp.data") {
		if (missing(snpsubset)) snpsubset=data@snpnames
		if (missing(idsubset))  idsubset=data@idnames
		n1data <- data[idsubset,snpsubset]
	} else {
		stop("data argument should be of type gwaa.data or snp.data");
	}
	if (missing(hweidsubset)) hweidsubset <- n1data@idnames

	nts <- n1data@nsnps
	cat(nts,"markers in total\n")
	out <- list()
# check redundancy
#	browser()
	redok <- rep(1,n1data@nsnps)
	if (redundant != "no" && (redundant == "bychrom" || redundant=="all")) {
		out$details.redundancy <- redundant(n1data,pairs=redundant, minconcordance = minconcordance)
		out$redundant <- out$details.redundancy[["all"]]
		redok[match(out$details.redundancy$all,n1data@snpnames)] <- 0
	}
	ntmp <- length(out$redundant); ptmp <- 100*ntmp/nts
	cat(ntmp," (",ptmp,"%) markers excluded as redundant\n",sep="")
# run summary
	s <- summary(n1data)
# check frequency
	mrmaf <- pmin(s$Q.2,1-s$Q.2)
	freqok <- (mrmaf>=maf)
	out$nofreq <- n1data@snpnames[which(freqok == FALSE)]
	ntmp <- length(out$nofreq); ptmp <- 100*ntmp/nts
	cat(ntmp," (",ptmp,"%) markers excluded as having low maf\n",sep="")
# check call rate
	callok <- (s[,"CallRate"]>=callrate)
	out$nocall <- n1data@snpnames[which(callok == FALSE)]
	ntmp <- length(out$nocall); ptmp <- 100*ntmp/nts
	cat(ntmp," (",ptmp,"%) markers excluded because of low call rate\n",sep="")
# check HWE
	shwe <- summary(n1data[hweidsubset,])
	Pv <- shwe$Pexact
	if (p.level < 0) {
		if (qoption == "storey") {
			qobj <- qvalue(Pv,fdr.level=fdrate)
			hweok <- !(qobj$significant)
		} else {
			hweok <- !(qvaluebh95(Pv,fdrate)$significant)
		}
	} else {
		hweok <- (Pv >= p.level)
	}
	out$nohwe <- n1data@snpnames[which(hweok == FALSE)]
#	cat("nohwe: ",out$nohwe,"\n")
	ntmp <- length(out$nohwe); ptmp <- 100*ntmp/nts
	cat(ntmp," (",ptmp,"%) markers excluded because they are out of HWE\n",sep="")
# all together
	out$ok <- n1data@snpnames[(callok & freqok & redok & hweok)]
#	cat("ok: ",out$ok,"\n")
	ntmp <- length(out$ok); ptmp <- 100*ntmp/nts
	cat("In total, ",ntmp," (",ptmp,"%) markers passed all criteria\n",sep="")
#Chi2 for the ones out of HWE, sorted
	out$Pex.nohwe <- s$Pexact[match(out$nohwe,n1data@snpnames)]
#	cat("Pex.nohwe: ",out$Pex.nohwe,"\n")
	idx <- sort(out$Pex.nohwe,dec=FALSE,ind=TRUE)$ix
	out$nohwe <- out$nohwe[idx]
#	cat("nohwe: ",out$nohwe,"\n")
	out$Pex.nohwe <- out$Pex.nohwe[idx]
#	cat("Pex.nohwe: ",out$Pex.nohwe,"\n")
# what was the call
	out$call$call <- match.call()
	out$call$name <- n1data@snpnames
	out$call$map <- n1data@map
	out$call$chromosome <- n1data@chromosome
# output
	class(out) <- "check.marker"
	out
}

