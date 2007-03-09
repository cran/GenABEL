"check.marker" <-
function(data, snpsubset, idsubset,
			callrate=0.95,perid.call=0.95,
			het.fdr=0.01, ibs.threshold = 0.95, ibs.mrk = 2000, maf, p.level=-1, 
			fdrate = 0.2, odds = 100, hweidsubset, redundant="no", minconcordance = 2.0, 
			qoption="bh95") {

	if (class(data) == "gwaa.data") data <- data@gtdata
	if (class(data) != "snp.data") stop("data argument should be of type gwaa.data or snp.data");
	
	Xidfail <- NULL
	Xmrkfail <- NULL
	if (any(data@chromosome=="X") & any(data@male==1)) {
		cat(data@nsnps,"markers and",data@nids,"people in total\n")
		cat("Running sex (X-chromosome) checks...\n")
		flush.console()
		xmrk <- (data@chromosome == "X")
		mlst <- (data@male == 1)
		rxm <- data[mlst,xmrk]
		Xout <- Xcheck(rxm)
		if (Xout$xerr) {
			outerr <- Xout$tab
			cat("Wrong male X genotypes (heterozygous) found for",dim(Xout$tab)[1],"genotypes\n")
			cat("Error table is saved in Xerrtab\n")
			flush.console()
			xerrsnps <- unique(outerr[,2])
			Pssw <- 0.01
			for (jj in xerrsnps) {
				P2 <- summary(data[,jj])[,"Q.2"]
				PnX <- as.vector((-1.)*as.numeric(data[mlst,jj]))
				PssX <- PnX
				notna <- !is.na(PnX)
				vec0 <- (PnX==0 & notna)
				vec1 <- (PnX==-1 & notna)
				vec2 <- (PnX==-2 & notna)
				PnX[vec0] <- (1.-P2)*(1.-P2) 
				PnX[vec1] <- 2.*P2*(1.-P2) 
				PnX[vec2] <- P2*P2
				PnonX <- sum(log(PnX),na.rm=TRUE)
				PssX[vec0] <- (1.-P2)*(1.-Pssw)
				PssX[vec1] <- 2.*(1.-P2)*P2*Pssw
				PssX[vec2] <- P2*(1.-Pssw)
				PseswX <- sum(log(PssX),na.rm=TRUE)
				if ((PnonX-PseswX)>log(1000)) {
					cat("Marker",jj,"is likely to be not an X marker (Odds>",odds,")\n")
					flush.console()
					if (is.null(Xmrkfail)) {
						Xmrkfail <- jj
					} else {
						Xmrkfail <- c(Xmrkfail,jj)
					}
				}
			}
			xerrids <- unique(outerr[,1])
			Pssw <- 0.01
			Pgterr <- 0.01
			for (jj in xerrids) {
				P2 <- summary(data[,xmrk])[,"Q.2"]
				PssX <- as.vector((-1.)*as.numeric(data[jj,xmrk]))
				PerrX <- PssX
				notna <- !is.na(PssX)
				vec0 <- (PssX==0 & notna)
				vec1 <- (PssX==-1 & notna)
				vec2 <- (PssX==-2 & notna)
				PssX[vec0] <- (1.-P2[vec0])*(1.-P2[vec0])
				PssX[vec1] <- 2.*P2[vec1]*(1.-P2[vec1])
				PssX[vec2] <- P2[vec2]*P2[vec2]
				PseswX <- sum(log(PssX),na.rm=TRUE) #+log(Pssw)
				PerrX[vec0] <- (1.-P2[vec0])*(1.-Pgterr)
				PerrX[vec1] <- Pgterr
				PerrX[vec2] <- P2[vec2]*(1.-Pgterr)
				PerrorX <- sum(log(PerrX),na.rm=TRUE)
				if ((PseswX-PerrorX)>log(1000)) {
					cat("Person",jj,"is likely to be female (Odds>",odds,")\n")
					flush.console()
					if (is.null(Xidfail)) {
						Xidfail <- jj
					} else {
						Xidfail <- c(Xidfail,jj)
					}
				}
			}
		} else {
			cat("No sex errors found\n")
			flush.console()
			outerr <- NULL
		}
		if (!is.null(Xidfail) | !is.null(Xmrkfail)) {
			if (!is.null(Xidfail)) Xidok <- rxm@idnames[!(rxm@idnames %in% Xidfail)] else Xidok <- rxm@idnames
			if (!is.null(Xmrkfail)) Xsnpok <- rxm@snpnames[!(rxm@snpnames %in% Xmrkfail)] else Xsnpok <- rxm@snpnames
			Xout <- Xcheck(rxm[Xidok,Xsnpok])
			cat("If these people / snps are removed ");
			if (Xout$xerr) {
				cat("wrong male X genotypes (heterozygous) found for",dim(Xout$tab)[1],"genotypes\n")
			} else {cat("no errors are detected\n")}
		}
		rm(rxm,xmrk,mlst);gc(verbose=FALSE)
	} else outerr <- NULL;
#	if (!is.null(Xidfail)) Xidok <- data@idnames[!(data@idnames %in% Xidfail)]
#	if (!is.null(Xmrkfail)) Xsnpok <- data@snpnames[!(data@snpnames %in% Xmrkfail)]
	run <- 1
	cat("RUN",run,"\n");
	if (!missing(hweidsubset)) {
		if (is.logical(hweidsubset)) hweidsubset <- data@idnames[which(hweidsubset==TRUE)]
		if (is.numeric(hweidsubset)) hweidsubset <- data@idnames[hweidsubset[!is.na(hweidsubset)]]
	}
	out0 <- check.marker.internal(data=data, snpsubset=snpsubset, idsubset=idsubset,
			callrate=callrate,perid.call=perid.call,
			het.fdr=het.fdr, ibs.threshold = ibs.threshold, ibs.mrk = ibs.mrk, maf = maf, p.level= p.level, 
			fdrate = fdrate, hweidsubset = hweidsubset, redundant = redundant, minconcordance = minconcordance, 
			qoption = qoption)
	out0$Xidfail <- Xidfail
	out0$Xmrkfail <- Xmrkfail
	out0$idok <- out0$idok[!(out0$idok %in% Xidfail)]
	out0$snpok <- out0$snpok[!(out0$snpok %in% Xmrkfail)]
	out1 <- list(); out1$idok <- ""; out1$snpok <- "";
	call <- out0$call
	run <- run + 1
	if (length(out0$idok)==data@nids && length(out0$snpok)==data@nsnps) {
		passed <- TRUE
		out1 <- out0 
	} else passed <- FALSE
	while (!passed) {
		cat("RUN",run,"\n");
		if (!missing(hweidsubset)) hweidsubset <- hweidsubset[hweidsubset %in% out0$idok]
		out1 <- check.marker.internal(data=data[out0$idok,out0$snpok], 
			callrate=callrate,perid.call=perid.call,
			het.fdr=het.fdr, ibs.threshold = ibs.threshold, ibs.mrk = ibs.mrk, maf = maf, p.level= p.level, 
			fdrate = fdrate, hweidsubset=hweidsubset, redundant = redundant, minconcordance = minconcordance, 
			qoption = qoption)
		if (length(out1$idok) == length(out0$idok))
		if (all(out1$idok == out0$idok))
		if (length(out1$snpok) == length(out0$snpok))
		if (all(out1$snpok == out0$snpok)) {
			passed <- TRUE
			out1$nohwe <- out0$nohwe
			out1$Pex.nohwe <- out0$Pex.nohwe
			out1$nocall <- out0$nocall
			out1$nofreq <- out0$nofreq
			out1$redundant <- out0$redundant
			out1$idnocall <- out0$idnocall
			out1$hetfail <- out0$hetfail
			out1$ibsfail <- out0$ibsfail
			out1$Xidfail <- out0$Xidfal
			out1$Xmrkfail <- out0$Xmrkfail
			out1$details.redundancy[["all"]] <- out0$details.redundancy[["all"]]
			for (i in names(out0$details.redundancy)) 
				if (i != "all") 
					out1$details.redundancy[[i]] <- out0$details.redundancy[[i]]
		}
		if (!passed) {
			out1$nohwe <- c(out0$nohwe,out1$nohwe)
			out1$Pex.nohwe <- c(out0$Pex.nohwe,out1$Pex.nohwe)
			out1$nocall <- c(out0$nocall,out1$nocall)
			out1$nofreq <- c(out0$nofreq,out1$nofreq)
			out1$redundant <- c(out0$redundant,out1$redundant)
			out1$idnocall <- c(out0$idnocall,out1$idnocall)
			out1$hetfail <- c(out0$hetfail,out1$hetfail)
			out1$ibsfail <- c(out0$ibsfail,out1$ibsfail)
			out1$Xidfail <- out0$Xidfal
			out1$Xmrkfail <- out0$Xmrkfail
			out1$details.redundancy[["all"]] <- c(out0$details.redundancy[["all"]],out1$details.redundancy[["all"]])
			for (i in names(out0$details.redundancy)) 
				if (i != "all") 
					out1$details.redundancy[[i]] <- out0$details.redundancy[[i]]
			out0 <- out1
			run <- run + 1
		}
	}
	out1$call <- call
	out1$Xidfail <- Xidfail
	out1$Xmrkfail <- Xmrkfail
	out1$Xerrtab <- outerr
	out1
}
