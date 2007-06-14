"Xfix" <-
function(data) {
	if (class(data) != "gwaa.data") stop("data argument must be of gwaa.data-class")
	if(!any(data@gtdata@chromosome=="X")) stop("no X-linked markers in the data")
	mlst <- data@gtdata@male==1
	xmrk <- data@gtdata@snpnames[data@gtdata@chromosome=="X"]
	err <- Xcheck(data@gtdata[mlst,xmrk])
	if (!err$xerr) {
		cat("no X-errors to fix\n")
		return(data)
	} else {
		for (i in 1:dim(err$tab)[1]) {
			cid <- which(data@gtdata@idnames==err$tab[i,1])
			csn <- which(data@gtdata@snpnames==err$tab[i,2])
			gtv <- as.numeric(data@gtdata[,csn])+1
			gtv[is.na(gtv)] <- 0
			gtv[cid] <- 0
			data@gtdata@gtps[,csn] <- put.snps(gtv)
		}
	}
	data
}
