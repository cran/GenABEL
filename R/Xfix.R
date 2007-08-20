"Xfix" <-
function(data) {
	if (class(data) != "gwaa.data") stop("data argument must be of gwaa.data-class")
	if(!any(data@gtdata@chromosome=="X")) stop("no X-linked markers in the data")
	mlst <- data@gtdata@male==1
	xmrk <- data@gtdata@chromosome=="X"
	if (sum(mlst)<1 || sum(xmrk)<1) {
		cat("no SNPs typed in male\n")
		return(data)
	}
	err <- Xcheck(data@gtdata[mlst,xmrk],tabonly=T)
	if (!err$xerr) {
		cat("no X-errors to fix\n")
		return(data)
	} else {
		for (i in 1:dim(err$Xerrtab)[1]) {
			cid <- which(data@gtdata@idnames==err$Xerrtab[i,1])
			csn <- which(data@gtdata@snpnames==err$Xerrtab[i,2])
			gtv <- as.numeric(data@gtdata[,csn])+1
			gtv[is.na(gtv)] <- 0
			gtv[cid] <- 0
			data@gtdata@gtps[,csn] <- put.snps(gtv)
		}
		cat(dim(err$Xerrtab)[1],"X-errors fixed\n")
	}
	data
}
