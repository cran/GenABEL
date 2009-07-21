"perid.summary" <- 
function (data,snpsubset,idsubset, ...) {
	if (class(data)=="gwaa.data") data <- data@gtdata
	if (class(data)!="snp.data") stop("The data argument must be of snp.data-class or gwaa.data-class")
	if (!missing(snpsubset)) data <- data[,snpsubset]
	if (!missing(idsubset)) data <- data[idsubset,]
	out <- hom(data, ... )
	out <- as.data.frame(out,stringsAsFactors=F)
	out$CallPP <- out$NoMeasured/data@nsnps
	out$Het <- 1. - out$Hom
	rownames(out) <- data@idnames
	out
}
