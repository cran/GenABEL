"ibs" <- 
function (data,snpsubset,idsubset,weight="no") {
	if (class(data)=="gwaa.data") data <- data@gtdata
	if (class(data)!="snp.data") stop("The data argument must be of snp.data-class or gwaa.data-class")
	if (!missing(snpsubset)) data <- data[,snpsubset]
	if (!missing(idsubset)) data <- data[idsubset,]
	wargs <- c("no","freq")
	if (!(match(weight,wargs,nomatch=0)>0)) {
		out <- paste("weight argument should be one of",wargs,"\n")
		stop(out)
	}
#	if (npairs > 2500000) stop("Too many pairs to evaluate... Stopped")
	if (weight == "no") option = 0
	if (weight == "freq") option = 1
	out <- .C("ibs", as.raw(data@gtps), as.integer(data@nids), as.integer(data@nsnps), as.integer(option), sout = double(data@nids*data@nids), PACKAGE="GenABEL")$sout
	dim(out) <- c(data@nids,data@nids)
	out <- t(out)
	if (weight=="no") {
		diag(out) <- hom(data,weight="no")[,2]
	} else {
		diag(out) <- 0.5+hom(data,weight="freq")[,4]
	}
	colnames(out) <- data@idnames
	rownames(out) <- data@idnames
	out
}
